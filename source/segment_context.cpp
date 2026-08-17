#include "segment_context.h"

#include "support_functions.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>

namespace radiseq {
namespace {

int weighted_index_1_based(
    const std::vector<double>& weights,
    std::mt19937& generator)
{
    if (weights.empty()) {
        return 0;
    }

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const double selected = distribution(generator);
    double cumulative = 0.0;

    for (std::size_t i = 0; i < weights.size(); ++i) {
        cumulative += weights[i];
        if (selected <= cumulative) {
            return static_cast<int>(i + 1);
        }
    }

    // Protect against a tiny floating-point normalization deficit.
    return static_cast<int>(weights.size());
}

template <typename Integer>
Integer uniform_integer(
    Integer lower,
    Integer upper,
    std::mt19937& generator)
{
    if (upper < lower) {
        throw std::invalid_argument("Invalid integer sampling interval");
    }

    std::uniform_int_distribution<Integer> distribution(lower, upper);
    return distribution(generator);
}

double uniform_real(
    double lower,
    double upper,
    std::mt19937& generator)
{
    if (upper < lower) {
        throw std::invalid_argument("Invalid real sampling interval");
    }

    std::uniform_real_distribution<double> distribution(lower, upper);
    return distribution(generator);
}

bool build_gc_bias(
    SegmentContext& context,
    int minimum_fragment_length)
{
    if (!context.gc_bias_enabled) {
        return true;
    }

    const std::size_t sequence_length = context.sequence.size();
    if (sequence_length < static_cast<std::size_t>(context.read_length)) {
        return false;
    }

    std::bitset<256> is_n;
    std::bitset<256> is_gc;
    is_n.set(static_cast<unsigned char>('N'));
    is_gc.set(static_cast<unsigned char>('G'));
    is_gc.set(static_cast<unsigned char>('C'));

    int gc_count = 0;
    int n_count = 0;
    int number_of_bins = 0;
    double total_gc_bias = 0.0;
    std::vector<int> forbidden_bin_ids;

    for (std::size_t i = 0; i < sequence_length; ++i) {
        const unsigned char base =
            static_cast<unsigned char>(context.sequence[i]);

        if (is_gc.test(base)) {
            ++gc_count;
        } else if (is_n.test(base)) {
            ++n_count;
        }

        if ((i + 1U) % static_cast<std::size_t>(context.gc_bin_size) == 0U) {
            double bin_bias = 0.0;
            const int non_n_bases = context.gc_bin_size - n_count;
            if (non_n_bases != 0) {
                bin_bias = GCBias::get_GCbias(
                    static_cast<double>(gc_count) /
                    static_cast<double>(non_n_bases));
            }

            context.gc_bins_bias.push_back(bin_bias);
            total_gc_bias += bin_bias;
            ++number_of_bins;
            gc_count = 0;
            n_count = 0;

            if (bin_bias == 0.0) {
                // Existing ART code stores one-based bin identifiers.
                forbidden_bin_ids.push_back(number_of_bins);
            }
        }
    }

    const int unprocessed_segment =
        static_cast<int>(sequence_length) -
        number_of_bins * context.gc_bin_size;

    if (unprocessed_segment >= context.read_length) {
        double bin_bias = 0.0;
        const int non_n_bases = unprocessed_segment - n_count;
        if (non_n_bases != 0) {
            bin_bias = GCBias::get_GCbias(
                static_cast<double>(gc_count) /
                static_cast<double>(non_n_bases));
        }

        context.gc_bins_bias.push_back(bin_bias);
        total_gc_bias += bin_bias;
        ++number_of_bins;

        if (bin_bias == 0.0) {
            forbidden_bin_ids.push_back(number_of_bins);
        }
    }

    if (total_gc_bias == 0.0 || context.gc_bins_bias.empty()) {
        return false;
    }

    for (int bin_id : forbidden_bin_ids) {
        context.forbidden_bins.insert(bin_id);
    }

    double cumulative_gc_bias = 0.0;
    context.cumulative_gc_bins_bias.reserve(context.gc_bins_bias.size());
    for (double& bin_bias : context.gc_bins_bias) {
        bin_bias /= total_gc_bias;
        cumulative_gc_bias += bin_bias;
        context.cumulative_gc_bins_bias.push_back(cumulative_gc_bias);
    }

    if (minimum_fragment_length >= static_cast<int>(sequence_length)) {
        if (context.forbidden_bins.find(1) != context.forbidden_bins.end()) {
            return false;
        }
    } else {
        const int available_bins =
            static_cast<int>(std::floor(
                static_cast<double>(
                    static_cast<int>(sequence_length) -
                    minimum_fragment_length) /
                static_cast<double>(context.gc_bin_size))) +
            1;

        if (available_bins <= 0 ||
            static_cast<std::size_t>(available_bins) >
                context.cumulative_gc_bins_bias.size() ||
            context.cumulative_gc_bins_bias[
                static_cast<std::size_t>(available_bins - 1)] == 0.0) {
            return false;
        }
    }

    return static_cast<int>(context.forbidden_bins.size()) != number_of_bins;
}

bool build_mda_distribution(
    SegmentContext& context,
    std::uint64_t requested_pairs,
    double coverage_per_cell,
    std::mt19937& generator)
{
    if (!context.mda_enabled) {
        return true;
    }

    if (requested_pairs == 0 || coverage_per_cell == 0.0) {
        return false;
    }

    std::uint64_t initial_site_count = requested_pairs;
    if (coverage_per_cell > 1.0) {
        initial_site_count = static_cast<std::uint64_t>(
            static_cast<double>(requested_pairs) / coverage_per_cell);
    }

    if (initial_site_count == 0 ||
        initial_site_count >
            static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        return false;
    }

    const std::size_t site_count = static_cast<std::size_t>(initial_site_count);
    context.initial_mda_sites.clear();
    context.initial_mda_sites.reserve(site_count);
    context.initial_mda_site_bias.assign(site_count, 0.0);

    for (std::size_t i = 0; i < site_count; ++i) {
        long initial_site = 0;

        if (!context.gc_bias_enabled) {
            initial_site = uniform_integer<long>(
                0L,
                context.valid_region,
                generator);
        } else {
            const int initial_bin =
                weighted_index_1_based(context.gc_bins_bias, generator);
            if (initial_bin <= 0) {
                return false;
            }

            const long bin_start =
                static_cast<long>(initial_bin - 1) * context.gc_bin_size;
            const long bin_end = std::min(
                context.valid_region,
                static_cast<long>(initial_bin) * context.gc_bin_size);

            initial_site = uniform_integer<long>(
                bin_start,
                bin_end,
                generator);
        }

        context.initial_mda_sites.push_back(initial_site);
    }

    std::sort(
        context.initial_mda_sites.begin(),
        context.initial_mda_sites.end());

    const long number_of_biased_sites = uniform_integer<long>(
        1L,
        static_cast<long>(context.initial_mda_sites.size()),
        generator);

    double maximum_bias = 0.4;
    double cumulative_probability = 0.0;

    for (long j = 0; j < number_of_biased_sites; ++j) {
        const std::size_t site_index = static_cast<std::size_t>(
            uniform_integer<long>(
                1L,
                static_cast<long>(site_count),
                generator) -
            1L);

        const double site_bias =
            uniform_real(0.0, maximum_bias, generator);

        context.initial_mda_site_bias[site_index] += site_bias;
        cumulative_probability += site_bias;
        maximum_bias -= site_bias;

        if (maximum_bias <= 0.0) {
            break;
        }
    }

    const double remaining_probability = 1.0 - cumulative_probability;
    double total_bias = 0.0;

    for (double& bias : context.initial_mda_site_bias) {
        if (remaining_probability != 0.0) {
            bias += remaining_probability /
                static_cast<double>(site_count);
        }
        total_bias += bias;
    }

    if (total_bias != 0.0) {
        for (double& bias : context.initial_mda_site_bias) {
            bias /= total_bias;
        }
    }

    return !context.initial_mda_site_bias.empty();
}

} // namespace

SegmentContextPtr build_segment_context(
    std::string chromosome_id,
    std::string chromosome_sequence,
    const SegmentBuildOptions& options)
{
    if (options.read_length <= 0) {
        throw std::invalid_argument("SegmentContext read length must be positive");
    }
    if (options.gc_bin_size <= 0) {
        throw std::invalid_argument("SegmentContext GC bin size must be positive");
    }
    if (options.minimum_fragment_length <= 0) {
        throw std::invalid_argument(
            "SegmentContext minimum fragment length must be positive");
    }

    auto context = std::make_shared<SegmentContext>();
    context->chromosome_id = std::move(chromosome_id);
    context->sequence = std::move(chromosome_sequence);
    context->read_length = options.read_length;
    context->gc_bin_size = options.gc_bin_size;
    context->valid_region =
        static_cast<long>(context->sequence.size()) - options.read_length;
    context->gc_bias_enabled = GCBias::get_GCbias_slope() != 0.0;
    context->mda_enabled = options.use_mda;

    if (context->valid_region < 1) {
        return {};
    }

    if (!build_gc_bias(*context, options.minimum_fragment_length)) {
        return {};
    }

    std::seed_seq seed_sequence{
        static_cast<unsigned int>(options.random_seed),
        static_cast<unsigned int>(options.random_seed >> 32U),
        static_cast<unsigned int>(options.segment_index),
        static_cast<unsigned int>(options.segment_index >> 32U),
        static_cast<unsigned int>(context->sequence.size()),
        static_cast<unsigned int>(context->sequence.size() >> 32U)
    };
    std::mt19937 loader_generator(seed_sequence);

    if (!build_mda_distribution(
            *context,
            options.requested_pairs,
            options.coverage_per_cell,
            loader_generator)) {
        return {};
    }

    return context;
}

} // namespace radiseq
