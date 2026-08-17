#ifndef RADISEQ_PAIRED_END_PIPELINE_H
#define RADISEQ_PAIRED_END_PIPELINE_H

#include "parameter_handler.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace radiseq {

struct PairedPipelineConfig {
    // Number of successful read pairs requested in each work item.
    std::size_t pairs_per_job = 4096;

    // Number of fully preprocessed chromosome contexts allowed ahead of the dispatcher.
    // Keep this small because each context owns an entire chromosome sequence.
    std::size_t segment_queue_capacity = 2;

    // Bounded queue capacities. A zero value selects a worker-count-based default.
    std::size_t work_queue_capacity = 0;
    std::size_t output_queue_capacity = 0;

    // zlib compression level: 1 is fastest, 6 is zlib's default, 9 is smallest.
    int gzip_compression_level = 3;

    // Existing RadiSeq code adds 32 to quality values. Keep 32 initially so
    // this refactor does not silently change output semantics. Change to 33
    // only after separately validating the quality-profile convention.
    int quality_ascii_offset = 32;

    // Abort a work item if repeated ART rejection makes it effectively impossible.
    std::size_t max_attempts_per_pair = 10000;

    // Preserve deterministic block order in the FASTQ files even though jobs
    // finish out of order.
    bool preserve_block_order = true;
};

struct PairedPipelineStats {
    std::uint64_t segments_loaded = 0;
    std::uint64_t segments_accepted = 0;
    std::uint64_t jobs_submitted = 0;
    std::uint64_t pairs_generated = 0;
    std::uint64_t generation_attempts = 0;
    std::uint64_t output_blocks_written = 0;
};

// Runs paired-end sequencing for one memory-mapped single-cell FASTA file.
// The loader builds immutable SegmentContext objects and may preprocess several
// chromosomes ahead while a dispatcher feeds chunked jobs to ART workers.
PairedPipelineStats run_single_cell_paired_pipeline(
    NGSParameters& parameter,
    const char* fasta_data,
    std::size_t fasta_size,
    double coverage_per_cell,
    int min_fragment_length,
    const std::vector<double>& fragment_weights,
    const std::string& output_r1_path,
    const std::string& output_r2_path,
    const PairedPipelineConfig& config = PairedPipelineConfig{});

} // namespace radiseq

#endif
