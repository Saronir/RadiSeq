#include "paired_end_pipeline.h"

#include "art_framework.h"
#include "bounded_queue.h"
#include "fileio.h"
#include "segment_context.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>
#include <zlib.h>

namespace radiseq {
namespace {

struct SegmentTask {
    SegmentContextPtr context;
    std::uint64_t requested_pairs = 0;
};

struct ReadJob {
    std::uint64_t block_id = 0;
    std::uint64_t first_pair_id = 0;
    std::size_t pair_count = 0;
    SegmentContextPtr segment;
};

struct FastqBlock {
    std::uint64_t block_id = 0;
    std::uint64_t pair_count = 0;
    std::string r1;
    std::string r2;
};

struct SharedFailure {
    void capture(std::exception_ptr error)
    {
        bool expected = false;
        if (failed.compare_exchange_strong(expected, true, std::memory_order_acq_rel)) {
            std::lock_guard<std::mutex> lock(mutex);
            first_error = std::move(error);
        }
    }

    void rethrow_if_set() const
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (first_error) {
            std::rethrow_exception(first_error);
        }
    }

    std::atomic<bool> failed{false};
    mutable std::mutex mutex;
    std::exception_ptr first_error;
};

class PairedGzipWriter {
public:
    PairedGzipWriter(
        const std::string& r1_path,
        const std::string& r2_path,
        int compression_level)
    {
        r1_ = gzopen(r1_path.c_str(), "wb");
        r2_ = gzopen(r2_path.c_str(), "wb");

        if (r1_ == nullptr || r2_ == nullptr) {
            close_noexcept();
            throw std::runtime_error("Unable to open paired FASTQ gzip output files");
        }

        // Increase zlib's internal buffers to reduce small system calls.
        gzbuffer(r1_, 1U << 20);
        gzbuffer(r2_, 1U << 20);

        const int level = std::clamp(compression_level, 0, 9);
        if (gzsetparams(r1_, level, Z_DEFAULT_STRATEGY) != Z_OK ||
            gzsetparams(r2_, level, Z_DEFAULT_STRATEGY) != Z_OK) {
            close_noexcept();
            throw std::runtime_error("Unable to configure gzip compression");
        }
    }

    ~PairedGzipWriter()
    {
        close_noexcept();
    }

    PairedGzipWriter(const PairedGzipWriter&) = delete;
    PairedGzipWriter& operator=(const PairedGzipWriter&) = delete;

    void write(const FastqBlock& block)
    {
        write_all(r1_, block.r1);
        write_all(r2_, block.r2);
    }

private:
    static void write_all(gzFile file, std::string_view data)
    {
        while (!data.empty()) {
            const std::size_t max_chunk =
                static_cast<std::size_t>(std::numeric_limits<unsigned int>::max());
            const unsigned int chunk = static_cast<unsigned int>(
                std::min<std::size_t>(data.size(), max_chunk));

            const int written = gzwrite(file, data.data(), chunk);
            if (written <= 0) {
                int error_code = Z_OK;
                const char* message = gzerror(file, &error_code);
                throw std::runtime_error(
                    message != nullptr ? message : "Unknown gzip write failure");
            }

            data.remove_prefix(static_cast<std::size_t>(written));
        }
    }

    void close_noexcept() noexcept
    {
        if (r1_ != nullptr) {
            gzclose(r1_);
            r1_ = nullptr;
        }
        if (r2_ != nullptr) {
            gzclose(r2_);
            r2_ = nullptr;
        }
    }

    gzFile r1_ = nullptr;
    gzFile r2_ = nullptr;
};

std::string clean_chromosome_id(std::string id)
{
    if (!id.empty() && id.front() == '>') {
        id.erase(id.begin());
    }
    return id;
}

void append_quality(
    std::string& output,
    const std::vector<short>& quality,
    std::size_t sequence_length,
    int ascii_offset)
{
    if (quality.size() < sequence_length) {
        throw std::runtime_error("ART returned fewer quality values than read bases");
    }

    for (std::size_t i = 0; i < sequence_length; ++i) {
        const int encoded = static_cast<int>(quality[i]) + ascii_offset;
        if (encoded < 0 || encoded > 126) {
            throw std::runtime_error("FASTQ quality value is outside the printable ASCII range");
        }
        output.push_back(static_cast<char>(encoded));
    }
}

void append_fastq_record(
    std::string& output,
    const std::string& chromosome_id,
    std::uint64_t pair_id,
    int mate_number,
    const std::string& sequence,
    const std::vector<short>& quality,
    int quality_ascii_offset)
{
    output.push_back('@');
    output.append(chromosome_id);
    output.append("_read");
    output.append(std::to_string(pair_id));
    output.push_back('/');
    output.push_back(static_cast<char>('0' + mate_number));
    output.push_back('\n');
    output.append(sequence);
    output.append("\n+\n");
    append_quality(output, quality, sequence.size(), quality_ascii_offset);
    output.push_back('\n');
}

} // namespace

PairedPipelineStats run_single_cell_paired_pipeline(
    NGSParameters& parameter,
    const char* fasta_data,
    std::size_t fasta_size,
    double coverage_per_cell,
    int min_fragment_length,
    const std::vector<double>& fragment_weights,
    const std::string& output_r1_path,
    const std::string& output_r2_path,
    const PairedPipelineConfig& config)
{
    if (fasta_data == nullptr || fasta_size == 0) {
        throw std::invalid_argument("Paired pipeline received an empty FASTA memory map");
    }
    if (config.pairs_per_job == 0) {
        throw std::invalid_argument("pairs_per_job must be greater than zero");
    }
    if (config.segment_queue_capacity == 0) {
        throw std::invalid_argument("segment_queue_capacity must be greater than zero");
    }
    if (config.max_attempts_per_pair == 0) {
        throw std::invalid_argument("max_attempts_per_pair must be greater than zero");
    }

    const std::size_t worker_count = static_cast<std::size_t>(
        std::max(1, parameter.get_number_of_threads()));
    const std::size_t work_capacity = config.work_queue_capacity != 0
        ? config.work_queue_capacity
        : std::max<std::size_t>(4, worker_count * 3);
    const std::size_t output_capacity = config.output_queue_capacity != 0
        ? config.output_queue_capacity
        : std::max<std::size_t>(4, worker_count * 2);

    BoundedQueue<SegmentTask> segment_queue(config.segment_queue_capacity);
    BoundedQueue<ReadJob> work_queue(work_capacity);
    BoundedQueue<FastqBlock> output_queue(output_capacity);
    SharedFailure failure;

    std::atomic<std::uint64_t> segments_loaded{0};
    std::atomic<std::uint64_t> segments_accepted{0};
    std::atomic<std::uint64_t> jobs_submitted{0};
    std::atomic<std::uint64_t> pairs_generated{0};
    std::atomic<std::uint64_t> generation_attempts{0};
    std::atomic<std::uint64_t> blocks_written{0};
    std::atomic<std::size_t> workers_remaining{worker_count};

    ART read1;
    ART read2;
    read1.initialize_worker_storage(static_cast<int>(worker_count));
    read2.initialize_worker_storage(static_cast<int>(worker_count));
    read1.set_read_error_rates(
        parameter.get_insertion_error_rate_read1(),
        parameter.get_deletion_error_rate_read1());
    read2.set_read_error_rates(
        parameter.get_insertion_error_rate_read2(),
        parameter.get_deletion_error_rate_read2());

    auto close_all_queues = [&] {
        segment_queue.close();
        work_queue.close();
        output_queue.close();
    };

    std::thread writer([&] {
        try {
            PairedGzipWriter gzip_writer(
                output_r1_path,
                output_r2_path,
                config.gzip_compression_level);

            std::uint64_t next_block = 0;
            std::map<std::uint64_t, FastqBlock> pending;
            FastqBlock block;

            while (output_queue.pop(block)) {
                if (!config.preserve_block_order) {
                    gzip_writer.write(block);
                    blocks_written.fetch_add(1, std::memory_order_relaxed);
                    continue;
                }

                pending.emplace(block.block_id, std::move(block));
                auto it = pending.find(next_block);
                while (it != pending.end()) {
                    gzip_writer.write(it->second);
                    pending.erase(it);
                    blocks_written.fetch_add(1, std::memory_order_relaxed);
                    ++next_block;
                    it = pending.find(next_block);
                }
            }

            if (config.preserve_block_order &&
                !pending.empty() &&
                !failure.failed.load(std::memory_order_acquire)) {
                throw std::runtime_error(
                    "Writer finished with a missing FASTQ block");
            }
        } catch (...) {
            failure.capture(std::current_exception());
            close_all_queues();
        }
    });

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (std::size_t worker_id = 0; worker_id < worker_count; ++worker_id) {
        workers.emplace_back([&, worker_id] {
            try {
                std::vector<short> quality1;
                std::vector<short> quality2;
                quality1.reserve(
                    static_cast<std::size_t>(parameter.get_read_length()));
                quality2.reserve(
                    static_cast<std::size_t>(parameter.get_read_length()));

                ReadJob job;
                while (!failure.failed.load(std::memory_order_acquire) &&
                       work_queue.pop(job)) {
                    if (!job.segment) {
                        throw std::runtime_error(
                            "Worker received a null SegmentContext");
                    }

                    FastqBlock block;
                    block.block_id = job.block_id;
                    block.pair_count = job.pair_count;

                    const std::size_t estimated_record_bytes =
                        static_cast<std::size_t>(parameter.get_read_length()) *
                            2U +
                        96U;
                    block.r1.reserve(job.pair_count * estimated_record_bytes);
                    block.r2.reserve(job.pair_count * estimated_record_bytes);

                    std::size_t successful = 0;
                    std::size_t attempts = 0;
                    const std::size_t max_attempts =
                        job.pair_count * config.max_attempts_per_pair;

                    while (successful < job.pair_count) {
                        if (failure.failed.load(std::memory_order_acquire)) {
                            break;
                        }
                        if (++attempts > max_attempts) {
                            throw std::runtime_error(
                                "ART rejected too many paired-read generation attempts for " +
                                job.segment->chromosome_id);
                        }

                        generation_attempts.fetch_add(
                            1,
                            std::memory_order_relaxed);

                        if (!ART::generate_paired_reads_with_indel(
                                *job.segment,
                                read1,
                                read2,
                                min_fragment_length,
                                fragment_weights,
                                static_cast<int>(worker_id))) {
                            continue;
                        }

                        read1.get_read_quality(
                            quality1,
                            1,
                            static_cast<int>(worker_id));
                        read1.add_baseCall_error(
                            quality1,
                            static_cast<int>(worker_id));
                        const std::string& sequence1 =
                            *read1.get_final_read_sequence(
                                static_cast<int>(worker_id));

                        read2.get_read_quality(
                            quality2,
                            2,
                            static_cast<int>(worker_id));
                        read2.add_baseCall_error(
                            quality2,
                            static_cast<int>(worker_id));
                        const std::string& sequence2 =
                            *read2.get_final_read_sequence(
                                static_cast<int>(worker_id));

                        const std::uint64_t pair_id =
                            job.first_pair_id + successful;
                        append_fastq_record(
                            block.r1,
                            job.segment->chromosome_id,
                            pair_id,
                            1,
                            sequence1,
                            quality1,
                            config.quality_ascii_offset);
                        append_fastq_record(
                            block.r2,
                            job.segment->chromosome_id,
                            pair_id,
                            2,
                            sequence2,
                            quality2,
                            config.quality_ascii_offset);

                        ++successful;
                    }

                    if (successful != job.pair_count) {
                        break;
                    }

                    pairs_generated.fetch_add(
                        successful,
                        std::memory_order_relaxed);
                    if (!output_queue.push(std::move(block))) {
                        break;
                    }
                }
            } catch (...) {
                failure.capture(std::current_exception());
                close_all_queues();
            }

            if (workers_remaining.fetch_sub(
                    1,
                    std::memory_order_acq_rel) == 1) {
                output_queue.close();
            }
        });
    }

    std::thread dispatcher([&] {
        try {
            std::uint64_t next_block_id = 0;
            std::uint64_t next_pair_id = 0;
            SegmentTask segment_task;

            while (!failure.failed.load(std::memory_order_acquire) &&
                   segment_queue.pop(segment_task)) {
                if (!segment_task.context) {
                    throw std::runtime_error(
                        "Dispatcher received a null SegmentContext");
                }

                std::uint64_t submitted_pairs = 0;
                while (submitted_pairs < segment_task.requested_pairs) {
                    const std::uint64_t remaining =
                        segment_task.requested_pairs - submitted_pairs;
                    const std::size_t count = static_cast<std::size_t>(
                        std::min<std::uint64_t>(
                            remaining,
                            config.pairs_per_job));

                    ReadJob job;
                    job.block_id = next_block_id++;
                    job.first_pair_id = next_pair_id;
                    job.pair_count = count;
                    job.segment = segment_task.context;

                    if (!work_queue.push(std::move(job))) {
                        break;
                    }

                    jobs_submitted.fetch_add(1, std::memory_order_relaxed);
                    submitted_pairs += count;
                    next_pair_id += count;
                }
            }
        } catch (...) {
            failure.capture(std::current_exception());
            close_all_queues();
        }
        work_queue.close();
    });

    std::thread loader([&] {
        try {
            std::size_t position = 0;
            std::string chromosome_sequence;
            std::string chromosome_id;
            std::uint64_t segment_index = 0;

            while (!failure.failed.load(std::memory_order_acquire) &&
                   getNextChromSeq_MM(
                       fasta_data,
                       fasta_size,
                       position,
                       chromosome_sequence,
                       chromosome_id)) {
                segments_loaded.fetch_add(1, std::memory_order_relaxed);

                const long requested_reads = static_cast<long>(
                    (coverage_per_cell *
                     static_cast<double>(chromosome_sequence.size())) /
                    static_cast<double>(parameter.get_read_length()));

                if (requested_reads < 1 ||
                    chromosome_sequence.size() <=
                        static_cast<std::size_t>(
                            parameter.get_read_length())) {
                    ++segment_index;
                    continue;
                }

                const std::uint64_t requested_pairs =
                    (static_cast<std::uint64_t>(requested_reads) + 1U) /
                    2U;

                SegmentBuildOptions build_options;
                build_options.read_length = parameter.get_read_length();
                build_options.gc_bin_size = parameter.get_GC_binSize();
                build_options.minimum_fragment_length =
                    min_fragment_length;
                build_options.use_mda =
                    *parameter.get_coverage_distribution() == "mda";
                build_options.requested_pairs = requested_pairs;
                build_options.coverage_per_cell = coverage_per_cell;
                build_options.random_seed = parameter.get_random_seed();
                build_options.segment_index = segment_index++;

                SegmentContextPtr context = build_segment_context(
                    clean_chromosome_id(chromosome_id),
                    std::move(chromosome_sequence),
                    build_options);

                if (!context) {
                    continue;
                }

                segments_accepted.fetch_add(1, std::memory_order_relaxed);
                SegmentTask task;
                task.context = std::move(context);
                task.requested_pairs = requested_pairs;

                if (!segment_queue.push(std::move(task))) {
                    break;
                }
            }
        } catch (...) {
            failure.capture(std::current_exception());
            close_all_queues();
        }
        segment_queue.close();
    });

    loader.join();
    dispatcher.join();
    for (std::thread& worker : workers) {
        worker.join();
    }
    writer.join();

    failure.rethrow_if_set();

    PairedPipelineStats stats;
    stats.segments_loaded = segments_loaded.load();
    stats.segments_accepted = segments_accepted.load();
    stats.jobs_submitted = jobs_submitted.load();
    stats.pairs_generated = pairs_generated.load();
    stats.generation_attempts = generation_attempts.load();
    stats.output_blocks_written = blocks_written.load();
    return stats;
}

} // namespace radiseq
