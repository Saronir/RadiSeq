#ifndef FASTAFILE_HANDLER_H
#define FASTAFILE_HANDLER_H

#include <string>

#include "sddfile_handler.h"

//long buildUndamagedGenomeTemplate(const std::string&, int, int, const std::string*);
long buildUndamagedGenomeTemplate_MM(char*, std::size_t, int, int, const std::string*, std::vector<double>&, double*, int);
double getReverseComplementarySeq(const std::string&, std::string&, int GC_binSize=0);                               // GC_binSize is optional
//int buildDamagedCellGenome(NGSsdd&, const std::string&, const std::string&);
std::vector<double> buildDamagedCellGenome_from_MM(NGSsdd&, const std::string&, const std::string&, char*, size_t, long, int);
std::vector<double> buildMutatedCellGenome_from_MM(const std::string&, const std::string&, const std::string&, char*, size_t, long, NGSParameters&);
std::string generateDNA(int);
std::string reverseComplement(const std::string&);
struct Mutation_Metadata{
    std::string mutation_type = "none"; //none, longdel, balinv, baltrans, indel
    double normalized_position = 0.0;
    int position = 0; //position on chromosome
    int length = 0; //length of mutation
    int pair = 0; //chromosome pair for translocations
    std::string inordel = "none"; //none, in, del
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

#endif