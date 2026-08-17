#ifndef RADISEQ_SEGMENT_CONTEXT_H
#define RADISEQ_SEGMENT_CONTEXT_H

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace radiseq {

// Immutable chromosome-specific input used by paired-end ART workers.
// A ReadJob stores a shared_ptr<const SegmentContext>, so the sequence and its
// GC/MDA distributions remain alive until the last job for the segment ends.
struct SegmentContext {
    std::string chromosome_id;
    std::string sequence;

    int read_length = 0;
    int gc_bin_size = 0;
    long valid_region = 0;

    bool gc_bias_enabled = false;
    bool mda_enabled = false;

    std::vector<double> gc_bins_bias;
    std::vector<double> cumulative_gc_bins_bias;
    std::unordered_set<long> forbidden_bins;

    std::vector<long> initial_mda_sites;
    std::vector<double> initial_mda_site_bias;
};

using SegmentContextPtr = std::shared_ptr<const SegmentContext>;

struct SegmentBuildOptions {
    int read_length = 0;
    int gc_bin_size = 0;
    int minimum_fragment_length = 1;

    bool use_mda = false;
    std::uint64_t requested_pairs = 0;
    double coverage_per_cell = 0.0;

    std::uint64_t random_seed = 0;
    std::uint64_t segment_index = 0;
};

// Returns nullptr when the segment cannot produce reads under the requested
// read length, GC-bias, fragment-length, or MDA constraints.
SegmentContextPtr build_segment_context(
    std::string chromosome_id,
    std::string chromosome_sequence,
    const SegmentBuildOptions& options);

} // namespace radiseq

#endif
