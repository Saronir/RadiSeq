#include "fileio.h"
#include "support_functions.h"
#include "random_generator.h"
#include "mutationfile_handler.h"
#include "json.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <sys/mman.h>
#include <omp.h>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <filesystem>

//READS JSON FILE CONTAINING CHROMOSOME DATA FOR EXACT LOCATIONS OF MUTATIONS FOLLOWING THE STRUCTURE USED 
//IN THE CDATA_MATRIX AND MDATA_MATRIX OBJECTS

//-LL

// using json = nlohmann::json;

//FUNCTION TO CHECK FOR MUTATION FILE IN WORKING DIRECTORY
bool check_mutation_file(const std::string& filename) {
    return std::filesystem::exists(filename);
}

//FUNCTION TO GENERATE MUTATIONS GIVEN MUTATIONAL PARAMETERS
std::vector<Chromosome_Metadata> mutation_generator(NGSParameters& parameters){

    int num_chrom = parameters.get_number_chromo();
    double SVfreq = parameters.get_structural_variation_frequency();
    double mean_num_long_del = parameters.get_proportion_long_deletions() * SVfreq;         //mean number of long deletions
    double mean_num_bal_inv = parameters.get_proportion_balanced_inversions() * SVfreq;
    double mean_num_bal_trans = parameters.get_proportion_balanced_translocations() * SVfreq;
    double mean_num_indel = parameters.get_proportion_indels() * SVfreq;
    int num_long_del = rng::poisson_sample(mean_num_long_del);                              //number of long deletions in the cell
    int num_bal_inv = rng::poisson_sample(mean_num_bal_inv);
    int num_bal_trans = rng::poisson_sample(mean_num_bal_trans);
    int num_indel = rng::poisson_sample(mean_num_indel);

    //sum the mutations and then distribute them to unique chromosomes, making sure you don't run out of chromosomes
    //some loss of generality in doing this, but otherwise there are troubling end cases and logical sorting to deal with
    int total_mutations = num_long_del + num_bal_inv + (2.0 * num_bal_trans) + num_indel;
    std::vector<int> mutated_chromosomes = rng::non_unique_random(1, num_chrom, total_mutations);

    // for (int kk = 0; kk < mutated_chromosomes.size(); kk++){
    //     std::cout << mutated_chromosomes[kk] << std::endl;
    // }

    //construct chromosome data mutation matrix from params 
    /*
        chromosome_number | mutation_metadata 
    
    */
    std::vector<Chromosome_Metadata> cdata_matrix; 
    std::vector<Mutation_Metadata> mdat_matrix;

    for (int i = 0; i < num_chrom; i++){
        int counter = 0;
        for (int j = 0; j < total_mutations; j++){
            if (mutated_chromosomes[j] == i + 1){
                counter += 1;
            }
        }
        Chromosome_Metadata chrom;
        chrom.chromosome_number = i + 1;
        chrom.num_mutations = counter;
        cdata_matrix.push_back(chrom);
    }

    //std::cout << "check2" << std::endl;
    // std::cin.get();

    //Populate mutation matrix with all required params besides location 
    int length_LD_min = parameters.get_min_long_deletion_length();          //fetch min length of long deletion
    int length_LD_max = parameters.get_max_long_deletion_length();          //fetch max length of long deletion
    for (int i = 0; i < num_long_del; i++){
        Mutation_Metadata mut;
        mut.mutation_type = "longdel";
        mut.length = rng::int_sample(length_LD_min, length_LD_max);     //pick a random length for each long deletion
        mut.normalized_position = rng::double_sample_0to1();
        // Chromosome_Metadata chrom;
        // chrom.chromosome_number = mutated_chromosomes[i];
        // chrom.mutation_metadata = mut;
        // cdata_matrix.push_back(chrom);
        mdat_matrix.push_back(mut);
    }
    //std::cout << "check2a" << std::endl;
    int length_BI_min = parameters.get_min_balanced_inversion_length();
    int length_BI_max = parameters.get_max_balanced_inversion_length();
    for (int i = 0; i < num_bal_inv; i++){
        Mutation_Metadata mut;
        mut.mutation_type = "balinv";
        mut.length = rng::int_sample(length_BI_min, length_BI_max);
        mut.normalized_position = rng::double_sample_0to1();
        // Chromosome_Metadata chrom;
        // chrom.chromosome_number = mutated_chromosomes[i + num_long_del];
        // chrom.mutation_metadata = mut;
        // cdata_matrix.push_back(chrom);
        mdat_matrix.push_back(mut);
    }
    //std::cout << "check2b" << std::endl;
    int length_BT_min = parameters.get_min_balanced_translocation_length();
    int length_BT_max = parameters.get_max_balanced_translocation_length();
    for (int i = 0; i < num_bal_trans; i++){
        Mutation_Metadata mut;
        mut.mutation_type = "baltrans";
        mut.length = rng::int_sample(length_BT_min, length_BT_max);
        //Chromosome_Metadata chrom;
        //chrom.chromosome_number = mutated_chromosomes[i + num_long_del + num_bal_inv];
        mut.pair = mutated_chromosomes[i + num_long_del + num_bal_inv + num_bal_trans];
        //chrom.mutation_metadata = mut;
        //cdata_matrix.push_back(chrom);
        mdat_matrix.push_back(mut);

        Mutation_Metadata mut_pair;
        mut_pair.mutation_type = "baltrans";
        mut_pair.length = rng::int_sample(length_BI_min, length_BI_max);
        //Chromosome_Metadata chrom_pair;
        //chrom_pair.chromosome_number = mut.pair;
        mut_pair.pair = mutated_chromosomes[i + num_long_del + num_bal_inv];
        //hrom_pair.mutation_metadata = mut_pair;
        //cdata_matrix.push_back(chrom_pair);
        mdat_matrix.push_back(mut_pair);
    }
    //std::cout << "check2c" << std::endl;
    int length_ID_min = parameters.get_min_indel_length();
    int length_ID_max = parameters.get_max_indel_length();
    for (int i = 0; i < num_indel; i++){
        Mutation_Metadata mut;
        mut.mutation_type = "indel";
        mut.length = rng::int_sample(500, 2000);
        mut.normalized_position = rng::double_sample_0to1();
        int coinFlip = rng::int_sample(0,1);
        std::cout << coinFlip << std::endl;
        if (coinFlip == 0){
            mut.inordel = "in";
        }
        else if (coinFlip == 1){
            mut.inordel = "del";
        }
        else{
            //throw some exception here
        }
        //Chromosome_Metadata chrom;
        //chrom.chromosome_number = mutated_chromosomes[i + num_long_del + num_bal_inv + 2 * num_bal_trans];
        //chrom.mutation_metadata = mut;
        //cdata_matrix.push_back(chrom);
        mdat_matrix.push_back(mut);
    }
    //std::cout << "check2d" << std::endl;
    //shuffle up the mutation metadata matrix
    rng::shuffleVector(mdat_matrix);

    //put the mutations into the mutated chromosomes
    for (int kk = 0; kk < mdat_matrix.size(); kk++){
        std::cout << mdat_matrix[kk].mutation_type << " " << mdat_matrix[kk].inordel  << " " << mdat_matrix[kk].length << " " << mdat_matrix[kk].normalized_position << std::endl;
    }
    //std::cin.get();
    int hold = 0;
    for (int i = 0; i < num_chrom; i++){
        int mut_used = 0;
        if (cdata_matrix[i].num_mutations == 0){
            //nothing
            //std::cout << i+1 << std::endl;
        }
        else{
            for (int j = 0; j < cdata_matrix[i].num_mutations; j++){
                cdata_matrix[i].mdata_matrix.push_back(mdat_matrix[j + hold]);
                //std::cout << i+1 << std::endl;
                mut_used += 1;
            }
            hold += mut_used;
        }
    }
    if (hold == total_mutations){
        //nothing
    }
    else{
        std::cout << hold << " is not " << total_mutations << std::endl;
        throw std::invalid_argument("Data skewed, mutations not used.");
        //throw exception here cause something wrong
    }

    return cdata_matrix;

}

//FUCTION TO INPUT AND CREATE CDATA_MATRIX FROM INPUT FILE
Cell_Metadata input_mutation_file(const std::string& filepath){
    std::ifstream file(filepath);

    json j;
    file >> j;

    Cell_Metadata cell = j.get<Cell_Metadata>();

    return cell;
}

//FUNCTION TO OUTPUT MUTATION FILE FROM CDATA_MATRIX
void output_mutation_file(Cell_Metadata& cell, const std::string& outfilepath){
    json j = cell;

    std::ofstream file(outfilepath);
    file << j.dump(4);
}

//USING ARGUMENT DEPENDANT LOOKUP FOR THE JSON I/O
void to_json(json& j, const Mutation_Metadata& data){
    j = json{
        {"mutation type", data.mutation_type},
        {"normalized position", data.normalized_position},
        {"position", data.position},
        {"length", data.length},
        {"pair", data.pair},
        {"insertion or deletion", data.inordel},
        {"base pairs", data.base_pairs}
    };
}

void to_json(json& j, const Chromosome_Metadata& data){
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

void to_json(json& j, const Cell_Metadata& data){
    j = json{
        {"cell id", data.cell_id},
        {"chromosome metadata", data.cdata_matrix}
    };
}

void from_json(const json& j, Mutation_Metadata& data){
    j.at("mutation type").get_to(data.mutation_type);
    j.at("normalized position").get_to(data.normalized_position);
    j.at("position").get_to(data.position);
    j.at("length").get_to(data.length);
    j.at("pair").get_to(data.pair);
    j.at("insertion or deletion").get_to(data.inordel);
    j.at("base pairs").get_to(data.base_pairs);
}

void from_json(const json& j, Chromosome_Metadata& data){
    j.at("chromosome number").get_to(data.chromosome_number);
    j.at("number of mutations").get_to(data.num_mutations);
    j.at("chromosome A ID").get_to(data.chromA_id);
    j.at("chromosome A sequence").get_to(data.chromA_seq);
    j.at("chromosome B ID").get_to(data.chromB_id);
    j.at("chromosome B sequence").get_to(data.chromB_seq);
    j.at("mutation metadata").get_to(data.mdata_matrix);
}

void from_json(const json& j, Cell_Metadata& data){
    j.at("cell id").get_to(data.cell_id);
    j.at("chromosome metadata").get_to(data.cdata_matrix);
}