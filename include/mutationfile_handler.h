#ifndef MUTATIONFILE_HANDLER_H
#define MUTATIONFILE_HANDLER_H

#include <string>
#include <vector>

#include "sddfile_handler.h"
#include "json.hpp"

using json = nlohmann::json;

struct Mutation_Metadata{
    std::string mutation_type = "none"; //none, longdel, balinv, baltrans, indel
    double normalized_position = 0.0;
    int position = 0; //position on chromosome
    int length = 0; //length of mutation
    int pair = 0; //chromosome pair for translocations
    std::string inordel = "none"; //none, in, del
    std::string base_pairs = "";
};

struct Chromosome_Metadata{
    int chromosome_number = 0;
    int num_mutations = 0;
    std::string chromA_id = "";
    std::string chromA_seq = "";
    std::string chromB_id = "";
    std::string chromB_seq = "";
    std::vector<Mutation_Metadata> mdata_matrix;
};

struct Cell_Metadata{
    std::string cell_id = "";
    std::vector<Chromosome_Metadata> cdata_matrix;
};

void to_json(json&, const Mutation_Metadata&);
void to_json(json&, const Chromosome_Metadata&);
void to_json(json&, const Cell_Metadata&);
void from_json(const json&, Mutation_Metadata&);
void from_json(const json&, Chromosome_Metadata&);
void from_json(const json&, Cell_Metadata&);
Cell_Metadata input_mutation_file(const std::string&);
std::vector<Chromosome_Metadata> mutation_generator(NGSParameters&);
void output_mutation_file(Cell_Metadata&, const std::string&);
bool check_mutation_file(const std::string&);

#endif