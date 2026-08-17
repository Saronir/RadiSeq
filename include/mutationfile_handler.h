#ifndef MUTATIONFILE_HANDLER_H
#define MUTATIONFILE_HANDLER_H

#include <cstdint>
#include <string>
#include <vector>

#include "sddfile_handler.h"
#include "json.hpp"

using json = nlohmann::json;

struct Mutation_Metadata {
    std::string mutation_type = "none";  // none, longdel, balinv, baltrans, delins, indel
    double normalized_position = 0.0;
    int position = 0;                    // 0-based coordinate on chromosome
    int length = 0;                      // mutation length in base pairs
    int pair = 0;                        // 1-based partner chromosome for paired SV records
    std::uint64_t event_id = 0;          // shared by records belonging to one SV event
    std::string inordel = "none";        // none, in, del
    std::string base_pairs;
};

struct Chromosome_Metadata {
    int chromosome_number = 0;
    int num_mutations = 0;
    std::string chromA_id;
    std::string chromA_seq;
    std::string chromB_id;
    std::string chromB_seq;
    std::vector<Mutation_Metadata> mdata_matrix;
};

struct Cell_Metadata {
    std::string cell_id;
    std::vector<Chromosome_Metadata> cdata_matrix;
};

void to_json(json&, const Mutation_Metadata&);
void to_json(json&, const Chromosome_Metadata&);
void to_json(json&, const Cell_Metadata&);

void from_json(const json&, Mutation_Metadata&);
void from_json(const json&, Chromosome_Metadata&);
void from_json(const json&, Cell_Metadata&);

Cell_Metadata input_mutation_file(const std::string& filepath);
std::vector<Chromosome_Metadata> mutation_generator(NGSParameters& parameters);
void output_mutation_file(Cell_Metadata& cell, const std::string& outfilepath);
bool check_mutation_file(const std::string& filename);

#endif
