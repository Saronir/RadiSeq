#include "art_framework.h"

#include "fastafile_handler.h"
#include "random_generator.h"
#include "segment_context.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

void ART::initialize_worker_storage(int number_of_workers)
{
    if (number_of_workers < 1) {
        throw std::invalid_argument("ART worker count must be positive");
    }

    indel_map_vec.clear();
    read_seq_vec.clear();
    indel_map_vec.resize(static_cast<std::size_t>(number_of_workers));
    read_seq_vec.resize(static_cast<std::size_t>(number_of_workers));
}

bool ART::generate_paired_reads_with_indel(
    const radiseq::SegmentContext& segment,
    ART& read_1,
    ART& read_2,
    int minFragSize,
    const std::vector<double>& fragment_weights,
    int threadID)
{
    const long chromSegment_length =
        static_cast<long>(segment.sequence.length());
    const int context_read_length = segment.read_length;
    const int context_gc_bin_size = segment.gc_bin_size;

    if (context_read_length != read_length) {
        throw std::runtime_error(
            "SegmentContext read length does not match ART configuration");
    }
    if (context_gc_bin_size <= 0 || chromSegment_length <= context_read_length) {
        return false;
    }

    long fragment_length = 0;
    const int MAX_ITERATIONS = 1000;
    int bad_sampling = 0;

    while (true) {
        if (minFragSize >= chromSegment_length) {
            fragment_length = std::max(
                chromSegment_length - 1,
                static_cast<long>(context_read_length));
            break;
        }

        const int weighted_fragment_index =
            rng::weighted_rand_int(fragment_weights, threadID);
        if (weighted_fragment_index <= 0) {
            return false;
        }

        const int minRandomFragSize =
            minFragSize + weighted_fragment_index - 1;
        fragment_length = std::max(
            minRandomFragSize,
            context_read_length);

        const long available_region =
            chromSegment_length - fragment_length;
        const long frag_bin = static_cast<long>(std::floor(
            static_cast<double>(available_region) /
            static_cast<double>(context_gc_bin_size)));

        if (segment.gc_bias_enabled) {
            if (frag_bin < 0 ||
                static_cast<std::size_t>(frag_bin) >=
                    segment.cumulative_gc_bins_bias.size()) {
                ++bad_sampling;
            } else if (segment.cumulative_gc_bins_bias[
                           static_cast<std::size_t>(frag_bin)] == 0.0) {
                ++bad_sampling;
            } else if (fragment_length < chromSegment_length) {
                break;
            } else {
                ++bad_sampling;
            }
        } else if (fragment_length < chromSegment_length) {
            break;
        } else {
            ++bad_sampling;
        }

        if (bad_sampling > MAX_ITERATIONS) {
            fragment_length = minFragSize;
            if (fragment_length >= chromSegment_length) {
                fragment_length = chromSegment_length - 1;
            }
            break;
        }
    }

    const int length_changed_1 = read_1.get_indel_map(threadID);
    const int length_changed_2 = read_2.get_indel_map(threadID);

    int adjusted_length_changed_1 = length_changed_1;
    int adjusted_length_changed_2 = length_changed_2;

    if ((context_read_length - adjusted_length_changed_1) > fragment_length) {
        adjusted_length_changed_1 =
            read_1.get_balanced_indel_map(threadID);
    }
    if ((context_read_length - adjusted_length_changed_2) > fragment_length) {
        adjusted_length_changed_2 =
            read_2.get_balanced_indel_map(threadID);
    }

    long start_position_1 = 0;
    long start_position_2 = 0;

    if (segment.mda_enabled) {
        if (segment.initial_mda_sites.empty() ||
            segment.initial_mda_site_bias.empty()) {
            return false;
        }

        const int amplicon_selection = rng::weighted_rand_int(
            segment.initial_mda_site_bias,
            threadID);
        if (amplicon_selection <= 0) {
            return false;
        }

        const std::size_t amplicon_index =
            static_cast<std::size_t>(amplicon_selection - 1);
        if (amplicon_index >= segment.initial_mda_sites.size()) {
            return false;
        }

        const long amplicon_start =
            segment.initial_mda_sites[amplicon_index];
        const long available_region =
            chromSegment_length - fragment_length;
        const long amplicon_end = std::min(
            available_region,
            amplicon_start + 12000L);

        if (amplicon_start > amplicon_end) {
            start_position_1 = static_cast<long>(std::floor(
                static_cast<double>(available_region) *
                rng::rand_double(0.0, 1.0, threadID)));
        } else {
            start_position_1 = rng::rand_int(
                amplicon_start,
                amplicon_end,
                threadID);
        }
    } else if (segment.gc_bias_enabled) {
        if (segment.gc_bins_bias.empty() ||
            segment.cumulative_gc_bins_bias.empty()) {
            return false;
        }

        const long available_region =
            chromSegment_length - fragment_length;
        const long num_available_fullBins = static_cast<long>(std::floor(
            static_cast<double>(available_region) /
            static_cast<double>(context_gc_bin_size)));
        const long remaining_chromLength =
            chromSegment_length -
            num_available_fullBins * context_gc_bin_size;
        const int bins_in_remaining_chrmLength =
            static_cast<int>(std::floor(
                static_cast<double>(remaining_chromLength) /
                static_cast<double>(context_gc_bin_size)));
        const long length_of_xtraBits =
            remaining_chromLength - fragment_length;

        double bias_in_xtraBits = 0.0;
        if (length_of_xtraBits >= context_read_length) {
            const std::size_t extra_bin_index =
                static_cast<std::size_t>(num_available_fullBins);
            if (extra_bin_index >= segment.gc_bins_bias.size()) {
                return false;
            }

            if (bins_in_remaining_chrmLength == 0) {
                if (remaining_chromLength <= 0) {
                    return false;
                }
                bias_in_xtraBits =
                    (segment.gc_bins_bias[extra_bin_index] /
                     static_cast<double>(remaining_chromLength)) *
                    static_cast<double>(length_of_xtraBits);
            } else {
                bias_in_xtraBits =
                    (segment.gc_bins_bias[extra_bin_index] /
                     static_cast<double>(context_gc_bin_size)) *
                    static_cast<double>(length_of_xtraBits);
            }
        }

        double max_cumGCbinBias = bias_in_xtraBits;
        if (num_available_fullBins != 0) {
            const std::size_t full_bin_index =
                static_cast<std::size_t>(num_available_fullBins - 1);
            if (full_bin_index >=
                segment.cumulative_gc_bins_bias.size()) {
                return false;
            }
            max_cumGCbinBias +=
                segment.cumulative_gc_bins_bias[full_bin_index];
        } else if (available_region == 0) {
            max_cumGCbinBias =
                segment.cumulative_gc_bins_bias.front();
        }

        if (max_cumGCbinBias <= 0.0) {
            return false;
        }

        const double random_binBias = rng::rand_double(
            0.0,
            max_cumGCbinBias,
            threadID);
        const auto selected = std::lower_bound(
            segment.cumulative_gc_bins_bias.begin(),
            segment.cumulative_gc_bins_bias.end(),
            random_binBias);

        if (selected == segment.cumulative_gc_bins_bias.end()) {
            return false;
        }

        const long selected_bin = static_cast<long>(std::distance(
            segment.cumulative_gc_bins_bias.begin(),
            selected));

        long binStart = 0;
        if (selected_bin <= num_available_fullBins) {
            binStart = selected_bin * context_gc_bin_size;
        } else {
            binStart = num_available_fullBins * context_gc_bin_size;
        }

        const long binEnd = std::min(
            available_region,
            (selected_bin + 1L) * context_gc_bin_size);
        if (binStart > binEnd) {
            return false;
        }

        start_position_1 = rng::rand_int(
            binStart,
            binEnd,
            threadID);
    } else {
        const long available_region =
            chromSegment_length - fragment_length;
        start_position_1 = static_cast<long>(std::floor(
            static_cast<double>(available_region) *
            rng::rand_double(0.0, 1.0, threadID)));
    }

    const int read_2_template_length = std::max(
        context_read_length,
        context_read_length - adjusted_length_changed_2);
    start_position_2 =
        (start_position_1 + fragment_length) - read_2_template_length;

    if (start_position_1 < 0 || start_position_2 < 0) {
        return false;
    }

    const long bin_1 =
        (start_position_1 / context_gc_bin_size) + 1;
    const long bin_2 =
        (start_position_2 / context_gc_bin_size) + 1;

    if (segment.forbidden_bins.find(bin_1) !=
            segment.forbidden_bins.end() ||
        segment.forbidden_bins.find(bin_2) !=
            segment.forbidden_bins.end()) {
        return false;
    }

    const int read_1_template_length = std::max(
        context_read_length,
        context_read_length - adjusted_length_changed_1);

    if (start_position_1 + read_1_template_length > chromSegment_length ||
        start_position_2 + read_2_template_length > chromSegment_length) {
        return false;
    }

    std::string read1_forward_seq = segment.sequence.substr(
        static_cast<std::size_t>(start_position_1),
        static_cast<std::size_t>(read_1_template_length));
    std::string read2_forward_seq = segment.sequence.substr(
        static_cast<std::size_t>(start_position_2),
        static_cast<std::size_t>(read_2_template_length));

    std::string read1_template_seq;
    std::string read2_template_seq;

    int read_orientation = 1;
    if (read_pair_orientation_dist[0] != 1.0) {
        read_orientation = rng::weighted_rand_int(
            read_pair_orientation_dist,
            threadID);
    }

    const int is_read1_chimeric = rng::weighted_rand_int(
        read_artifacts_rate_dist,
        threadID);
    const int is_read2_chimeric = rng::weighted_rand_int(
        read_artifacts_rate_dist,
        threadID);

    switch (read_orientation) {
        case 1:
            read1_template_seq = read1_forward_seq;
            getReverseComplementarySeq(
                read2_forward_seq,
                read2_template_seq);
            if (is_read1_chimeric == 1) {
                read_1.generate_chimeric_read(
                    read1_template_seq,
                    threadID);
            }
            if (is_read2_chimeric == 1) {
                read_2.generate_chimeric_read(
                    read2_template_seq,
                    threadID);
            }
            break;
        case 2:
            getReverseComplementarySeq(
                read1_forward_seq,
                read1_template_seq);
            read2_template_seq = read2_forward_seq;
            break;
        case 3:
            read1_template_seq = read1_forward_seq;
            read2_template_seq = read2_forward_seq;
            break;
        default:
            return false;
    }

    read_1.read_maker(read1_template_seq, threadID);
    read_2.read_maker(read2_template_seq, threadID);
    return true;
}
