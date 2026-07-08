#include "mutationfile_handler.h"

#include "random_generator.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void validate_proportion(double value, const std::string& name) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw std::invalid_argument(name + " must be between 0 and 1.");
    }
}

void validate_length_range(int minimum, int maximum, const std::string& name) {
    if (minimum <= 0 || maximum < minimum) {
        throw std::invalid_argument(
            name + " requires 0 < minimum <= maximum."
        );
    }
}

void add_mutation(
    std::vector<Chromosome_Metadata>& chromosomes,
    int chromosome_number,
    Mutation_Metadata mutation)
{
    if (chromosome_number < 1 ||
        chromosome_number > static_cast<int>(chromosomes.size()))
    {
        throw std::out_of_range("Mutation chromosome number is out of range.");
    }

    auto& chromosome = chromosomes.at(
        static_cast<std::size_t>(chromosome_number - 1)
    );

    chromosome.mdata_matrix.push_back(std::move(mutation));
    chromosome.num_mutations = static_cast<int>(chromosome.mdata_matrix.size());
}

} // namespace

bool check_mutation_file(const std::string& filename) {
    return !filename.empty() &&
           std::filesystem::is_regular_file(std::filesystem::path(filename));
}

std::vector<Chromosome_Metadata> mutation_generator(NGSParameters& parameters) {
    const int num_chrom = parameters.get_number_chromo();
    const double mutation_frequency =
        parameters.get_structural_variation_frequency();

    if (num_chrom <= 0) {
        throw std::invalid_argument("number_chromo must be positive.");
    }
    if (!std::isfinite(mutation_frequency) || mutation_frequency < 0.0) {
        throw std::invalid_argument(
            "structural_variation_frequency must be finite and non-negative."
        );
    }

    const double proportion_long_deletions =
        parameters.get_proportion_long_deletions();
    const double proportion_balanced_inversions =
        parameters.get_proportion_balanced_inversions();
    const double proportion_balanced_translocations =
        parameters.get_proportion_balanced_translocations();
    const double proportion_indels = parameters.get_proportion_indels();

    validate_proportion(proportion_long_deletions,
                        "proportion_long_deletions");
    validate_proportion(proportion_balanced_inversions,
                        "proportion_balanced_inversions");
    validate_proportion(proportion_balanced_translocations,
                        "proportion_balanced_translocations");
    validate_proportion(proportion_indels, "proportion_indels");

    const double proportion_sum =
        proportion_long_deletions +
        proportion_balanced_inversions +
        proportion_balanced_translocations +
        proportion_indels;

    if (proportion_sum > 1.0 + 1.0e-12) {
        throw std::invalid_argument(
            "The mutation-type proportions must sum to no more than 1."
        );
    }

    const int length_LD_min = parameters.get_min_long_deletion_length();
    const int length_LD_max = parameters.get_max_long_deletion_length();
    const int length_BI_min = parameters.get_min_balanced_inversion_length();
    const int length_BI_max = parameters.get_max_balanced_inversion_length();
    const int length_BT_min = parameters.get_min_balanced_translocation_length();
    const int length_BT_max = parameters.get_max_balanced_translocation_length();
    const int length_ID_min = parameters.get_min_indel_length();
    const int length_ID_max = parameters.get_max_indel_length();

    const int num_long_del = rng::poisson_sample(
        proportion_long_deletions * mutation_frequency
    );
    const int num_bal_inv = rng::poisson_sample(
        proportion_balanced_inversions * mutation_frequency
    );
    const int num_bal_trans = rng::poisson_sample(
        proportion_balanced_translocations * mutation_frequency
    );
    const int num_indel = rng::poisson_sample(
        proportion_indels * mutation_frequency
    );

    if (num_long_del > 0) {
        validate_length_range(
            length_LD_min, length_LD_max, "long deletion length"
        );
    }
    if (num_bal_inv > 0) {
        validate_length_range(
            length_BI_min, length_BI_max, "inversion length"
        );
    }
    if (num_bal_trans > 0) {
        validate_length_range(
            length_BT_min, length_BT_max, "translocation length"
        );
    }
    if (num_indel > 0) {
        validate_length_range(
            length_ID_min, length_ID_max, "indel length"
        );
    }

    if (num_bal_trans > 0 && num_chrom < 2) {
        throw std::invalid_argument(
            "At least two chromosomes are required for a translocation."
        );
    }

    std::vector<Chromosome_Metadata> chromosomes(
        static_cast<std::size_t>(num_chrom)
    );

    for (int i = 0; i < num_chrom; ++i) {
        chromosomes[static_cast<std::size_t>(i)].chromosome_number = i + 1;
    }

    for (int i = 0; i < num_long_del; ++i) {
        Mutation_Metadata mutation;
        mutation.mutation_type = "longdel";
        mutation.length = rng::int_sample(length_LD_min, length_LD_max);
        mutation.normalized_position = rng::double_sample_0to1();

        add_mutation(
            chromosomes,
            rng::int_sample(1, num_chrom),
            std::move(mutation)
        );
    }

    for (int i = 0; i < num_bal_inv; ++i) {
        Mutation_Metadata mutation;
        mutation.mutation_type = "balinv";
        mutation.length = rng::int_sample(length_BI_min, length_BI_max);
        mutation.normalized_position = rng::double_sample_0to1();

        add_mutation(
            chromosomes,
            rng::int_sample(1, num_chrom),
            std::move(mutation)
        );
    }

    std::uint64_t next_event_id = 1;

    for (int i = 0; i < num_bal_trans; ++i) {
        const int chromosome_a = rng::int_sample(1, num_chrom);

        int chromosome_b = chromosome_a;
        while (chromosome_b == chromosome_a) {
            chromosome_b = rng::int_sample(1, num_chrom);
        }

        const std::uint64_t event_id = next_event_id++;

        Mutation_Metadata mutation_a;
        mutation_a.mutation_type = "baltrans";
        mutation_a.length = rng::int_sample(length_BT_min, length_BT_max);
        mutation_a.pair = chromosome_b;
        mutation_a.event_id = event_id;

        Mutation_Metadata mutation_b;
        mutation_b.mutation_type = "baltrans";
        mutation_b.length = rng::int_sample(length_BT_min, length_BT_max);
        mutation_b.pair = chromosome_a;
        mutation_b.event_id = event_id;

        add_mutation(chromosomes, chromosome_a, std::move(mutation_a));
        add_mutation(chromosomes, chromosome_b, std::move(mutation_b));
    }

    for (int i = 0; i < num_indel; ++i) {
        Mutation_Metadata mutation;
        mutation.mutation_type = "indel";
        mutation.length = rng::int_sample(length_ID_min, length_ID_max);
        mutation.normalized_position = rng::double_sample_0to1();
        mutation.inordel = rng::int_sample(0, 1) == 0 ? "in" : "del";

        add_mutation(
            chromosomes,
            rng::int_sample(1, num_chrom),
            std::move(mutation)
        );
    }

    return chromosomes;
}

Cell_Metadata input_mutation_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Could not open mutation file: " + filepath);
    }

    json j;
    file >> j;
    return j.get<Cell_Metadata>();
}

void output_mutation_file(
    Cell_Metadata& cell,
    const std::string& outfilepath)
{
    std::ofstream file(outfilepath);
    if (!file) {
        throw std::runtime_error(
            "Could not open mutation output file: " + outfilepath
        );
    }

    const json j = cell;
    file << j.dump(4) << '\n';
}

void to_json(json& j, const Mutation_Metadata& data) {
    j = json{
        {"mutation type", data.mutation_type},
        {"normalized position", data.normalized_position},
        {"position", data.position},
        {"length", data.length},
        {"pair", data.pair},
        {"event id", data.event_id},
        {"insertion or deletion", data.inordel},
        {"base pairs", data.base_pairs}
    };
}

void to_json(json& j, const Chromosome_Metadata& data) {
    j = json{
        {"chromosome number", data.chromosome_number},
        {"number of mutations", data.num_mutations},
        {"chromosome A ID", data.chromA_id},
        {"chromosome A sequence", data.chromA_seq},
        {"chromosome B ID", data.chromB_id},
        {"chromosome B sequence", data.chromB_seq},
        {"mutation metadata", data.mdata_matrix}
    };
}

void to_json(json& j, const Cell_Metadata& data) {
    j = json{
        {"cell id", data.cell_id},
        {"chromosome metadata", data.cdata_matrix}
    };
}

void from_json(const json& j, Mutation_Metadata& data) {
    j.at("mutation type").get_to(data.mutation_type);
    data.normalized_position = j.value("normalized position", 0.0);
    data.position = j.value("position", 0);
    j.at("length").get_to(data.length);
    data.pair = j.value("pair", 0);
    data.event_id = j.value("event id", std::uint64_t{0});
    data.inordel = j.value("insertion or deletion", std::string{"none"});
    data.base_pairs = j.value("base pairs", std::string{});
}

void from_json(const json& j, Chromosome_Metadata& data) {
    j.at("chromosome number").get_to(data.chromosome_number);
    data.chromA_id = j.value("chromosome A ID", std::string{});
    data.chromA_seq = j.value("chromosome A sequence", std::string{});
    data.chromB_id = j.value("chromosome B ID", std::string{});
    data.chromB_seq = j.value("chromosome B sequence", std::string{});
    j.at("mutation metadata").get_to(data.mdata_matrix);
    data.num_mutations = static_cast<int>(data.mdata_matrix.size());
}

void from_json(const json& j, Cell_Metadata& data) {
    data.cell_id = j.value("cell id", std::string{});
    j.at("chromosome metadata").get_to(data.cdata_matrix);
}
