#include "mutationfile_handler.h"
#include "fastafile_handler.h"
#include "fileio.h"
#include "support_functions.h"
#include "random_generator.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sys/mman.h>
#include <omp.h>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include <filesystem>







//--------------------------------------------------------------------------------------------
// This function will read the reference genome file provided and generate an undamaged cell
// genome template that will be further used to create damaged cell genomes. The template will
// have forward and complementary strand sequences. If the cell is diploid, then it will even 
// have two copies of each strand. They will get IDs: copy1_chr1 and copy2_chr1 repectively for 
// two copies. Function also returns the total length of the reference sequence in bp.
//--------------------------------------------------------------------------------------------
/* long buildUndamagedGenomeTemplate(const std::string& filepath, int nChrms, int chrmMapping, const std::string* ref_seqPath){
    std::ofstream templateFile(filepath+"/Undamaged_cell.fa");                                            // File to store forward strand sequence. If you modify the name here, change in in other function below as well
    long seqLength{0};                                                                                    // Find the total length of the genome in bp from ref_seq file    
    if (templateFile.is_open()) {
        std::ifstream refseqfile(*ref_seqPath);
        int TotalChrms = nChrms;                                                                          // Holds the total number of chromosome till the end
        std::string chromSeq_ID;                                                                          // String to hold the chromosome sequence ID. This value is not used in this function but needed to pass to the getNextChromSeq functon
        std::string chromSeq;                                                                             // String to store the forward chromseq from reference seq file
        std::string revComp_chromSeq;                                                                     // String to store the reverse complementary sequence for respective chrom-seq
        int chrmCount{0};                                                                                 // Seperate counter to set the seq ID. This value will not be same as nChrms if diploid
        const int batchSize = 4000;                                                                       // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
        std::vector<std::string> batch_buffer;                                                            // Create a buffer for storing output data.
            
        switch(chrmMapping){                                                                              // Decide how to write the sequences to file depending on the chromosome mapping
            case 0:                                                                                       // If cell is haploid (chromosome mapping type 0)
                while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    getReverseComplementarySeq(chromSeq, revComp_chromSeq);                               // Generate the reverse complementary sequence of the chrom sequence
                    chrmCount++;
                    if(nChrms>2){                                                                         // Write according to the mapping 1 pattern
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--; 
                    }else{                                                                                // No need to have copies of X, Y chromosomes
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    }
                    if (batch_buffer.size() >= batchSize) {                                               // Check if the batch buffer is full, and write it to the file if needed.
                        writeBatchToFile(batch_buffer, templateFile);
                    }
                }
                break;
            case 1:                                                                                       // If the cell is diploid with chromosome mapping type 1 (1,1,2,2,....22,22,X,Y)
                while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    getReverseComplementarySeq(chromSeq, revComp_chromSeq);                               // Generate the reverse complementary sequence of the chrom sequence
                    chrmCount++;
                    if(nChrms>2){                                                                         // Until all autosomes are done, 
                        for(int i=0; i<2; i++){                                                           // Write according to the mapping 1 pattern
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                        } 
                    }else{                                                                                // No need to have copies of X, Y chromosomes
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    }
                    if (batch_buffer.size() >= batchSize) {                                               // Check if the batch buffer is full, and write it to the file if needed.
                        writeBatchToFile(batch_buffer, templateFile);
                    }
                }
                break;
            case 2:                                                                                       // If the cell is diploid with chromosome mapping type 2 (1,2.....22,1,2.....22,X,Y)
                for(int i=0; i<2; i++){
                    refseqfile.clear();                                                                   // Clear any error flags
                    refseqfile.seekg(0, std::ios::beg);                                                   // Set the position to the beginning of the file
                    chrmCount = 0; 
                    while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                        uppercaseString(chromSeq);                                                        // Change lowercase -> Uppercase
                        getReverseComplementarySeq(chromSeq, revComp_chromSeq);                           // Generate the reverse complementary sequence of the chrom sequence
                        chrmCount++;
                        if(chrmCount<int(TotalChrms/2)){                                                  // Write the autosomes once in the for loop
                            batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                            batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                            seqLength+=chromSeq.size();
                            nChrms--;
                        }
                        if (chrmCount>=int(TotalChrms/2) && i==0){                                        // Skip the sex chromosomes in the first for loop and write autosomes again 
                            break;
                        }else if(chrmCount>=int(TotalChrms/2) && i>0){                                    // In the second loop, write the sex chromosomes as well
                            batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                            batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                            seqLength+=chromSeq.size();
                            nChrms--;
                        }
                        if (batch_buffer.size() >= batchSize) {                                           // Check if the batch buffer is full, and write it to the file if needed.
                            writeBatchToFile(batch_buffer, templateFile);
                        }
                    }
                } 
                break;
        }
        if(nChrms!=0){                                                                                    // If fewer sequences than expected was found in the reference file, then error and exit
            std::cerr<<"\n ERROR: Only "<<chrmCount<<" chromosomse seqences were found in "<<*ref_seqPath<<"\n";
            std::cerr<<" A sequence file with "<<TotalChrms<<" chromosomes arranged in 1,2,3.....X,Y fashion is expected \n";
            exit(EXIT_FAILURE);
        }
        writeBatchToFile(batch_buffer, templateFile);                                                     // If there are unwritten data in batch buffer, write that too when the loop ends
        refseqfile.close();
    }
    templateFile.close();
    return seqLength;
} */
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will process the reference genome memory-map data provided and generate an 
// undamaged cell genome template that will be further used to create damaged cell genomes. 
// The template will have forward and complementary strand sequences. If the cell is diploid, 
// then it will even have two copies of each strand. They will get IDs: copy1_chr1 and copy2_chr1
// repectively for two copies. Function also returns the total length of the reference sequence
// in bp. Function will also populate a memory map of the file for easy handling.
//--------------------------------------------------------------------------------------------
long buildUndamagedGenomeTemplate_MM(char* templateFileMapping, std::size_t templateFileSize, int nChrms, int chrmMapping, const std::string* ref_seqPath, std::vector<double>& ref_chrm_weights, double* average_GC_content, int GC_binSize){
    char* position_in_MM = templateFileMapping;                                                           // Pointer to the current position in the MM as we write. Starts with the pointer to the beginning
    if (position_in_MM == nullptr){return 0;}                                                             // Return if the template file memory map is not valid
    
    long seqLength{0};                                                                                    // Find the total length of the genome in bp from ref_seq file  
    int TotalChrms = nChrms;                                                                              // Holds the total number of chromosome till the end
    std::string chromSeq_ID;                                                                              // String to hold the chromosome sequence ID. This value is not used in this function but needed to pass to the getNextChromSeq functon
    std::string chromSeq;                                                                                 // String to store the forward chromseq from reference seq file
    std::string revComp_chromSeq;                                                                         // String to store the reverse complementary sequence for respective chrom-seq
    int chrmCount{0};                                                                                     // Seperate counter to set the seq ID. This value will not be same as nChrms if diploid
    std::vector<long> chrm_length_weight;                                                                 // Temporary vector to hold the chromosome sizes before converting it to a weight
    std::vector<double> chrm_GC_fraction;                                                                 // Temporary vector to hold the ratio of GC count per chrom and the chrom size for each chromosome segment

    size_t position{0};                                                                                   // Temporary variable to hold the last-read position in the memory map of the reference file
    size_t refFileSize;                                                                                   // A variable to hold the file size of the reference genome, during memory-mapping
    void* refFileMM = generateInputFileMemoryMap(*ref_seqPath, refFileSize);                              // Create the memory-map of the reference genome file
    const char* refSeqData = static_cast<char*>(refFileMM);                                               // Casting the memory-map void pointer to a const char pointer for further processing
    const int batchSize{10};                                                                              // Define a batch size for writing data to the output memory-mapped file. These much data (buffer vector elements) will be stored in cache before writing it on the file
    std::vector<std::string> batch_buffer;                                                                // Create a buffer for storing output data.

    switch(chrmMapping){                                                                                  // Decide how to write the sequences to file depending on the chromosome mapping
        case 0:                                                                                           // If cell is haploid (chromosome mapping type 0)
            while(getNextChromSeq_MM(refSeqData, refFileSize, position, chromSeq, chromSeq_ID) && nChrms>0){
                uppercaseString(chromSeq);                                                                // Change lowercase -> Uppercase
                double GC_fraction = getReverseComplementarySeq(chromSeq, revComp_chromSeq, GC_binSize);   // Generate the reverse complementary sequence of the chrom sequence and also get the average GC fraction calculated for bins of given size
                chrmCount++;
                if(nChrms>2){                                                                             // Write according to the mapping 1 pattern
                    batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a\n"+chromSeq+"\n");
                    batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b\n"+revComp_chromSeq+"\n");
                    long chromSeqSize = static_cast<long>(chromSeq.size());
                    seqLength+=chromSeqSize;                                                              // Increment the reference seq length
                    chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);          // Append to the vector the genome size added in each fasta line. 0 corresponds to the chrom ID line in fasta
                    chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);
                    chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);             // Append to the vector the GC fraction for each line of fasta. 0 corresponds to the chrom ID line in fasta
                    chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);
                    nChrms--; 
                }else{                                                                                    // No need to have copies of X, Y chromosomes
                    batch_buffer.push_back(">chrXY_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                    batch_buffer.push_back(">chrXY_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                    long chromSeqSize = static_cast<long>(chromSeq.size());
                    seqLength+=chromSeqSize;                                                              // Increment the reference seq length
                    chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);          // Append to the vector the genome size added in each fasta line. 0 corresponds to the chrom ID line in fasta
                    chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);
                    chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);             // Append to the vector the GC fraction for each line of fasta. 0 corresponds to the chrom ID line in fasta
                    chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);
                    nChrms--;
                }
                if (batch_buffer.size() >= batchSize){                                                    // Check if the batch buffer is full, and write it to the memory-mapped file if needed.
                    writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);
                }    
            }
            break;
        case 1:                                                                                           // If the cell is diploid with chromosome mapping type 1 (1,1,2,2,....22,22,X,Y)
            while(getNextChromSeq_MM(refSeqData, refFileSize, position, chromSeq, chromSeq_ID) && nChrms>0){
                uppercaseString(chromSeq);                                                                // Change lowercase -> Uppercase
                double GC_fraction = getReverseComplementarySeq(chromSeq, revComp_chromSeq, GC_binSize);   // Generate the reverse complementary sequence of the chrom sequence and also get the average GC fraction calculated of read-sized bins
                chrmCount++;
                if(nChrms>2){                                                                             // Until all autosomes are done, 
                    for(int i=0; i<2; i++){                                                               // Write according to the mapping 1 pattern
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                        long chromSeqSize = static_cast<long>(chromSeq.size());
                        seqLength+=chromSeqSize;                                                          // Increment the reference seq length
                        chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);      // Append to the vector the genome size added in each fasta line. 0 corresponds to the chrom ID line in fasta
                        chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);
                        chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);         // Append to the vector the GC fraction for each line of fasta. 0 corresponds to the chrom ID line in fasta
                        chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);
                        nChrms--;
                    } 
                }else{                                                                                    // No need to have copies of X, Y chromosomes
                    batch_buffer.push_back(">chrXY_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                    batch_buffer.push_back(">chrXY_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                    long chromSeqSize = static_cast<long>(chromSeq.size());
                    seqLength+=chromSeqSize;                                                              // Increment the reference seq length
                    chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);          // Append to the vector the genome size added in each fasta line. 0 corresponds to the chrom ID line in fasta
                    chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);
                    chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);             // Append to the vector the GC fraction for each line of fasta. 0 corresponds to the chrom ID line in fasta
                    chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);
                    nChrms--;
                }
                if (batch_buffer.size() >= batchSize){                                                    // Check if the batch buffer is full, and write it to the memory-mapped file if needed.
                    writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);
                }  
            }
            break;
        case 2:                                                                                           // If the cell is diploid with chromosome mapping type 2 (1,2.....22,1,2.....22,X,Y)
            for(int i=0; i<2; i++){
                position = 0;                                                                             // Reset the starting positon of the memory map in each iteration
                chrmCount = 0; 
                while(getNextChromSeq_MM(refSeqData, refFileSize, position, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    double GC_fraction = getReverseComplementarySeq(chromSeq, revComp_chromSeq, GC_binSize);// Generate the reverse complementary sequence of the chrom sequence and also get the average GC fraction calculated of read-sized bins
                    chrmCount++;
                    if(chrmCount<int(TotalChrms/2)){                                                      // Write the autosomes once in the for loop
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                        long chromSeqSize = static_cast<long>(chromSeq.size());
                        seqLength+=chromSeqSize;                                                          // Increment the reference seq length
                        chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);      // Append to the vector the genome size added in each fasta line. 0 corresponds to the chrom ID line in fasta
                        chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);
                        chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);         // Append to the vector the GC fraction for each line of fasta. 0 corresponds to the chrom ID line in fasta
                        chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);
                        nChrms--;
                    }
                    if (chrmCount>=int((TotalChrms/2)-1) && i==0){                                        // Skip the sex chromosomes in the first for loop and write autosomes again 
                        break;
                    }else if(chrmCount>=int(TotalChrms/2) && i>0){                                        // In the second loop, write the sex chromosomes as well
                        batch_buffer.push_back(">chrXY_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chrXY_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                        long chromSeqSize = static_cast<long>(chromSeq.size());
                        seqLength+=chromSeqSize;                                                          // Increment the reference seq length
                        chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);      // Append to the vector the genome size added in each fasta line. 0 corresponds to the chrom ID line in fasta
                        chrm_length_weight.push_back(0); chrm_length_weight.push_back(chromSeqSize);
                        chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);         // Append to the vector the GC fraction for each line of fasta. 0 corresponds to the chrom ID line in fasta
                        chrm_GC_fraction.push_back(0.0); chrm_GC_fraction.push_back(GC_fraction);
                        nChrms--;
                    }
                    if (batch_buffer.size() >= batchSize){                                                // Check if the batch buffer is full, and write it to the memory-mapped file if needed.
                        writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);
                    }
                }
            } 
            break;
    }
    if(nChrms!=0){                                                                                        // If fewer sequences than expected was found in the reference file, then error and exit
        std::cerr<<"\n ERROR: Only "<<chrmCount<<" chromosomse seqences were found in the reference file\n";
        std::cerr<<" A sequence file with "<<TotalChrms<<" chromosomes arranged in 1,2,3.....X,Y fashion is expected \n";
        exit(EXIT_FAILURE);
    } 
    munmap(refFileMM, refFileSize);                                                                       // Unmap the reference genome file to avoid memory-leaks
    writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);              // If there are unwritten data in batch buffer, write that too when the loop ends                 
    
    *average_GC_content = averageOfEverySecond(chrm_GC_fraction);                                         // Calculate the average GC content as the average fraction of G,C in the entire genome
    GCBias::set_GCbias_peak(*average_GC_content);                                                         // Set where the GC bias to peak with its triangular function
    ref_chrm_weights.resize(chrm_length_weight.size());                                                   // Reshape the ref_chrm_weights vector to the required size in advance
    std::fill(ref_chrm_weights.begin(), ref_chrm_weights.end(), 0.0);                                     // Fill all the positions of the vector with zeros
    double total_chrm_weight{0.0};                                                                        // A temporary variable to hold the total weight from each chrm weight for normalization 
    #pragma omp parallel
    {
        #pragma omp for reduction (+:total_chrm_weight)
        for (size_t i = 1; i<chrm_length_weight.size(); i += 2){                                          // Skip the elements corresponding the chromosome IDs
            double length_weight = chrm_length_weight[i]/static_cast<double>(seqLength * 2);              // Divide each chromosome length with the total genome length from both strands combined to get the weights
            double GC_fraction_weight = GCBias::get_GCbias(chrm_GC_fraction[i]);                          // Get the GC bias from the triangular function according to the GC fraction
            ref_chrm_weights[i] = (length_weight*GC_fraction_weight);                                     // Populate the non-ID positions
            total_chrm_weight += (length_weight*GC_fraction_weight);                                      // Sum up all the individual chrm weights for normalization
        }

        #pragma omp for
        for (size_t i =1; i<ref_chrm_weights.size(); i += 2){                                             // Skip the elements corresponding to the chromosome IDs
            ref_chrm_weights[i] = ref_chrm_weights[i]/total_chrm_weight;                                  // Divide each chromosome length with the total weight to normalize the vector
        }
    }
    return seqLength;
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will make a reverse complementary sequence for the chromSeq referenced. The
// generated sequence will be in the direction 3'-5' if the original was 5'-3'. This new sequence
// will be stored in the reference string object called revComp_chromSeq. It will also count
// the number of G and C bases in the chromSeq passed and will return the mean GC fraction if 
// the third optional argument is specified. 
//--------------------------------------------------------------------------------------------
double getReverseComplementarySeq(const std::string& chromSeq, std::string& revComp_chromSeq, int GC_binSize){
    revComp_chromSeq.clear();                                                                             // Clear the string if there is previous seq
    revComp_chromSeq.resize(chromSeq.size());                                                             // Resize as with the forward sequence
    double GCfractions_allBins{0};                                                                        // Temporary variable to hold the sum of GC fractions in each read-sized bins
    long bin_counter{0};                                                                                  // Counter variable to count the number of bins with GC_binSize in the chromSeq passed
    int GC_counter{0};                                                                                    // Counter variable to count the number of G, C in the each bin with GC_binSize
    int N_counter{0};                                                                                     // Counter variable to count the number of Ns in each GC_binSize bin, to remove it from the bin size when calculating GC fraction

/*     std::ofstream outFile("GC_fractions.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return 0.0; // Return an error value
    } */
    for(size_t i=0; i<chromSeq.size(); i++){
        size_t k=chromSeq.size()-i-1;                                                                     // Indexing is reversed to get reverse strand
        switch(chromSeq[i]){
            case 'A':                                                                                     // Generate complementary base values accordingly
                revComp_chromSeq[k]='T'; break;
            case 'C':
                GC_counter++;
                revComp_chromSeq[k]='G'; break;
            case 'G':
                GC_counter++;
                revComp_chromSeq[k]='C'; break;
            case 'T':
                revComp_chromSeq[k]='A'; break;
            default:                                                                                      // Any unrecognized nitrogen base including 'N' will be set to 'N'
                N_counter++;
                revComp_chromSeq[k]='N';
        }

        if(GC_binSize != 0 && (i+1)%GC_binSize==0){                                                       // If i+1 index is a multiple of the GC_binSize. Bin size is zero by defualt. Skip the caculation if zero. 
            if (GC_binSize == N_counter){                                                                 // If all the elements in the bin are N's then GC fraction is zero
                GCfractions_allBins = 0.0;
            }else{
                GCfractions_allBins += (static_cast<double>(GC_counter)/(GC_binSize-N_counter));          // Add the GC fraction of the bin to the total. N count must be removed from GC_binSize for GC fraction calculation 
            }
            //           outFile <<(static_cast<double>(GC_counter)/read_size)<< '\t' <<GCBias::get_GCbias((static_cast<double>(GC_counter)/read_size))<< std::endl;
            bin_counter++; GC_counter = 0; N_counter = 0;                                                 // Increment the number of bins and reset the GC counter and N counter
        }
    }

    if (GC_binSize==0) return 0.0;                                                                        // No need to do the following calculation when GC binSize is zero
    int unprocessed_segment = chromSeq.size()-(bin_counter*GC_binSize);                                   // If there is more to the sequence after complete bins or if the sequence is smaller than one bin size
    if (unprocessed_segment>100){                                                                         // Check if this extra sequence is atleast greater than 100 bp, (100 is randomly chosen. Read length would be ideal)
        if (unprocessed_segment == N_counter){                                                            // If all the bases in the extra sequence are N's, GC fraction is zero
            GCfractions_allBins = 0.0;
        }else{
            GCfractions_allBins += (static_cast<double>(GC_counter)/(unprocessed_segment-N_counter));     // Add the GC fraction of the extra seq to the total. N count must be removed from GC_binSize for GC fraction calculation 
        }
        bin_counter++;                                                                                    // Increment the number of bins for this extra sequence for further calculations
    }
    return bin_counter!=0 ? (GCfractions_allBins/bin_counter):0.0;                                        // Retrun average GC fraction if bin_counter is not 0; else return 0
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will take the strand breaks, base damages and DNA breakpoint information obtained
// from SDD files and generate a damaged cell genome and store it into a fasta file. For every cell
// a fasta file will be made. The undamaged genome template previously generated will be used to 
// create the damaged cell genome. Base damage locations will be changed to an 'N' and DNA break
// points will be used to break the genome into segments. This function will return an integer
// value that is the total number of lines written in that particular fasta file for the damaged cell
//--------------------------------------------------------------------------------------------
/* int buildDamagedCellGenome(NGSsdd& SDDdata, const std::string& outputPath, const std::string& fileName) {
    std::ofstream outputFile(outputPath + fileName);
    if (!outputFile.is_open()) {return 0;}                                                                      // Return if the output file couldn't be opened

    std::ifstream genomeTemplate(outputPath + "/Undamaged_cell.fa");                                            // Undamaged genome template file that was previously built
    if (!genomeTemplate.is_open()) {return 0;}                                                                  // Return if the genome template file couldn't be opened

    int num_of_lines_written{0};                                                                                // This is a variable to hold the total number of lines written to each damage cell file
    std::string chromID_A;                                                                                      // Strings to hold the chrom IDs and Sequences from the file
    std::string chromSeq_A;
    std::string chromID_B;
    std::string chromSeq_B;
    long seqStartIndex = 0;                                                                                     // Variable that will hold the starting index for each chromsome sequence in the genome
    std::vector<long>& baseDamageLoc1 = SDDdata.get_basestrand1_damage_loc();                                   // Referencing the damage location vector to a temp vector for ease of handling
    std::vector<long>& baseDamageLoc2 = SDDdata.get_basestrand2_damage_loc();
    std::vector<long>& strandbreakLoc1 = SDDdata.get_backbone1_break_loc();
    std::vector<long>& strandbreakLoc2 = SDDdata.get_backbone2_break_loc();
    int i{0}; int j{0}; int k{0}; int l{0};                                                                     // Counter variables to keep a tab on the elements in the damage location vector that we already processed
    const int batchSize = 4000;                                                                                 // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
    std::vector<std::string> batch_buffer;                                                                      // Create a buffer for storing output data.
        
    while (readFastaTemplate(genomeTemplate, chromID_A, chromSeq_A, chromID_B, chromSeq_B)) {                   // Get forward, backward sequences and their repective IDs for each chroms, one at a time
        long seq_length = chromSeq_A.size();                                                                        // Sequence length is same for both A and B 

        // Introudce Base damages to the forward strand sequence
        while (i<baseDamageLoc1.size()) {                  
            if (baseDamageLoc1[i]<=(seqStartIndex + seq_length)) {                                              // If damage location is in the chromsome that is currently prcessing, then
                chromSeq_A[baseDamageLoc1[i] - seqStartIndex - 1] = 'N';                                        // Replace the nitrogen base at the damage location with 'N'. 1 is subtracted because the chrom_seq index starts from zero
                i++;                                                                                            // Increment the counter to know where the processing should start next in the damage location vector
            }else {break;}                                                                                      // Break the loop as soon as the damage location is more than the final base location of the chromosome. Works because the list is sorted. 
        }

        // Introudce Base damages to the reverse complementary strand sequence
        while (j<baseDamageLoc2.size()) {
            if (baseDamageLoc2[j]<=(seqStartIndex + seq_length)) {
                chromSeq_B[seq_length - (baseDamageLoc2[j] - seqStartIndex)] = 'N';                             // Indexing is different because the strand is not only complementary, but also revesrsed already. Index = chrom_end - (index in Genome - chrom_start)
                j++;
            } else {break;}
        }

        // Breaking chromosomes into segments, wherever there is a SSB
                
        // Processing the forward strand first 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        int A_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment in forward strand 
        if(strandbreakLoc1.size() !=0){                                                                         // If there are strand breaks present in the forward strand, then
            while(k<strandbreakLoc1.size()){
                if(strandbreakLoc1[k]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                    A_segment_length = (strandbreakLoc1[k]-seqStartIndex) - A_segment_start;
                    batch_buffer.push_back(chromSeq_A.substr(A_segment_start, A_segment_length)+"\n");                     // substr(a,b): get b charcters starting from index a 
                    A_segment_start = strandbreakLoc1[k]-seqStartIndex;
                    A_segment_index++; k++;
                    num_of_lines_written += 2;
                }else{                                                                                          // If the strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(A_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                        batch_buffer.push_back(chromSeq_A.substr(A_segment_start)+"\n");
                        num_of_lines_written += 2;
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_A+"\n");    
                        batch_buffer.push_back(chromSeq_A.substr(A_segment_start)+"\n");
                        num_of_lines_written += 2;
                    }
                    break;
                }
            }    
        }else{                                                                                                  // If there are no DNA breaks in the forward strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_A+"\n"+chromSeq_A+"\n");
            num_of_lines_written += 2;
        }
        
        // Processing the reverse complementary strand next 
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int B_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment of reverse strand
        if(strandbreakLoc2.size()!=0){                                                                          // If there are strand breaks present in reverse strand, then
            while(l<strandbreakLoc2.size()){                                                            
                if(strandbreakLoc2[l]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                    B_segment_length = B_segment_end - (seq_length-(strandbreakLoc2[l]-seqStartIndex));
                    batch_buffer.push_back(chromSeq_B.substr((seq_length-(strandbreakLoc2[l]-seqStartIndex)),B_segment_length)+"\n");
                    B_segment_end = seq_length-(strandbreakLoc2[l]-seqStartIndex);                            
                    B_segment_index++; l++;
                    num_of_lines_written += 2;
                }else{                                                                                          // If strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(B_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                        batch_buffer.push_back(chromSeq_B.substr(0,B_segment_end)+"\n");
                        num_of_lines_written += 2;
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_B+"\n");
                        batch_buffer.push_back(chromSeq_B.substr(0,B_segment_end)+"\n");
                        num_of_lines_written += 2;
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no strand breaks in the reverse strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_B+"\n"+chromSeq_B+"\n");
            num_of_lines_written += 2;
        }

        if (batch_buffer.size() >= batchSize) {                                                                 // Check if the batch buffer is full, and write it to the file if needed.
            writeBatchToFile(batch_buffer, outputFile);
        }
*/
        /* // Breaking chromosomes into segments, wherever there is a DSB 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int segment_index{1};                                                                                   // Temporary index variable to name each DNA segment  
        if(DNAbreaks.size()!=0){                                                                                // If there are DNA breaks present, then
            while(i<DNAbreaks.size()){                                                            
                if(DNAbreaks[i]<=(seqStartIndex+seq_length)){                                                   // If the DNA break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                    A_segment_length = (DNAbreaks[i]-seqStartIndex) - A_segment_start;
                    outputFile<<chromSeq_A.substr(A_segment_start, A_segment_length)<<"\n";                     // substr(a,b): get b charcters starting from index a 
                    A_segment_start = DNAbreaks[i]-seqStartIndex;

                    outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                    B_segment_length = B_segment_end - (seq_length-(DNAbreaks[i]-seqStartIndex));
                    outputFile<<chromSeq_B.substr((seq_length-(DNAbreaks[i]-seqStartIndex)),B_segment_length)<<"\n";
                    B_segment_end = seq_length-(DNAbreaks[i]-seqStartIndex);
                    
                    segment_index++;
                    DNAbreaks.erase(DNAbreaks.begin() + i);
                }else{                                                                                          // If DNA break is not in the chromosome Or if it is the final segment in a chromosome
                    if(segment_index!=1){                                                                       // if it is the last segment, write ID with segment tag
                        outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }else{                                                                                      // If the DNA break is not in chromosome, write without the segmenet tag
                        outputFile<<chromID_A<<"\n";    
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no DNA breaks in the genome, then just copy the sequence as is without segmenting
            outputFile<<chromID_A<<"\n"<<chromSeq_A<<"\n"<<chromID_B<<"\n"<<chromSeq_B<<"\n";
        } */
/*
        seqStartIndex += seq_length;
    }

    writeBatchToFile(batch_buffer, outputFile);                                                                 // If there are unwritten data in batch buffer, write that too when the loop ends

    genomeTemplate.close();
    outputFile.close();
    return (num_of_lines_written);
} */
//--------------------------------------------------------------------------------------------





//--------------------------------------------------------------------------------------------
// This function is almost identical to buildDamagedCellGenome function. The only difference is
// that this function uses the memory map of the undamaged genome template file to build the
// damaged genomes instead of using the file itself. If this function is called multiple times,
// this approach reduces the I/O operations on a large file hence improving performance. 
//--------------------------------------------------------------------------------------------
std::vector<double> buildDamagedCellGenome_from_MM(NGSsdd& SDDdata, const std::string& outputPath, const std::string& fileName, char* genomeTemplate_data, size_t templateSize, long ref_seq_length, int groupTID){
    size_t outFileSize = templateSize+10000;                                                                    // Temporary variable to hold the size of the current file being processed. This get dynamically changed if needed. Starts with undamagedfile size +10kb extra
    std::string outFilePath = outputPath+fileName;                                                              // Path to the output fasta file that we want to write
    char* outFileMapping = createMemoryMappedFile(outFilePath,outFileSize);                                     // Create a memory mapped output file for each damaged cell genome
    char* position_in_MM = outFileMapping;                                                                      // Pointer to the current position in the MM as we write. Starts with the pointer to the beginning
    std::vector<double> chrm_seg_weights;                                                                       // This is a vector to hold the weights of each chromosome segment scaled to its length
    if (position_in_MM == nullptr){return chrm_seg_weights;}                                                    // Return if the template file memory map is not valid
    std::vector<double> chrm_GC_bias;                                                                           // Temporary vector to hold the GC bias associated with each chromosome segment

    long total_seq_length = ref_seq_length*2;                                                                   // Twice the length because it has two strands contributing
    std::string chromID_A;                                                                                      // Strings to hold the chrom IDs and Sequences from the file
    std::string chromSeq_A;
    std::string chromID_B;
    std::string chromSeq_B;
    long seqStartIndex = 0;                                                                                     // Variable that will hold the starting index for each chromsome sequence in the genome
    std::vector<long>& baseDamageLoc1 = SDDdata.get_basestrand1_damage_loc(groupTID);                           // Referencing the damage location vector to a temp vector for ease of handling
    std::vector<long>& baseDamageLoc2 = SDDdata.get_basestrand2_damage_loc(groupTID);
    std::vector<long>& strandbreakLoc1 = SDDdata.get_backbone1_break_loc(groupTID);
    std::vector<long>& strandbreakLoc2 = SDDdata.get_backbone2_break_loc(groupTID);
    size_t i{0}; size_t j{0}; size_t k{0}; size_t l{0};                                                         // Counter variables to keep a tab on the elements in the damage location vector that we already processed
    const int batchSize = 100;                                                                                  // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
    std::vector<std::string> batch_buffer;                                                                      // Create a buffer for storing output data.
    size_t position = 0;                                                                                        // Temporary variable to hold the last read position in the memory map
    double GCslope = GCBias::get_GCbias_slope();                                                                // Temporary variable to hold the GC biase slope to avoid calling the function again and again
    

    while (readFastaMemoryMap(genomeTemplate_data, templateSize, position, chromID_A, chromSeq_A, chromID_B, chromSeq_B)){// Get forward, backward sequences and their repective IDs for each chroms, one at a time
        long seq_length = chromSeq_A.size();                                                                    // Sequence length is same for both A and B 
        
        // Introudce Base damages to the forward strand sequence
        while (i<baseDamageLoc1.size() && baseDamageLoc1[i] <= seqStartIndex + seq_length){                     // Proceed only if there are more unprocessed damages and if damage location is in the chromsome that is currently prcessing. Works because the list is sorted.
            chromSeq_A[baseDamageLoc1[i] - seqStartIndex - 1] = 'N';                                            // Replace the nitrogen base at the damage location with 'N'. 1 is subtracted because the chrom_seq index starts from zero
            i++;                                                                                                // Increment the counter to know where the processing should start next in the damage location vector
        }
        /* while (i<baseDamageLoc1.size()){                  
            if (baseDamageLoc1[i]<=(seqStartIndex + seq_length)){                                               // If damage location is in the chromsome that is currently prcessing, then
                chromSeq_A[baseDamageLoc1[i] - seqStartIndex - 1] = 'N';                                        // Replace the nitrogen base at the damage location with 'N'. 1 is subtracted because the chrom_seq index starts from zero
                i++;                                                                                            // Increment the counter to know where the processing should start next in the damage location vector
            }else {break;}                                                                                      // Break the loop as soon as the damage location is more than the final base location of the chromosome. Works because the list is sorted. 
        } */

        // Introudce Base damages to the reverse complementary strand sequence
        while (j < baseDamageLoc2.size() && baseDamageLoc2[j] <= seqStartIndex + seq_length){
            chromSeq_B[seq_length - (baseDamageLoc2[j] - seqStartIndex)] = 'N';                                 // Indexing is different because the strand is not only complementary, but also revesrsed already. Index = chrom_end - (index in Genome - chrom_start)
            ++j;
        }
        /* while (j<baseDamageLoc2.size()){
            if (baseDamageLoc2[j]<=(seqStartIndex + seq_length)){
                chromSeq_B[seq_length - (baseDamageLoc2[j] - seqStartIndex)] = 'N';                             // Indexing is different because the strand is not only complementary, but also revesrsed already. Index = chrom_end - (index in Genome - chrom_start)
                j++;
            } else {break;}
        } */

        // Breaking chromosomes into segments, wherever there is a SSB
                
        // Processing the forward strand first 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length{0};
        int A_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment in forward strand 
        if(k<strandbreakLoc1.size()){                                                                           // If there are unprocessed strand breaks present in the forward strand, then
            while(k<strandbreakLoc1.size()){
                if(strandbreakLoc1[k]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    A_segment_length = (strandbreakLoc1[k]-seqStartIndex) - A_segment_start;
                    if(A_segment_length==0){k++; continue;}                                                     // If the break is at the same location as the chromosome discontinuity, then continue to the next break
                    std::string chrm_segment_seq = chromSeq_A.substr(A_segment_start, A_segment_length);        // substr(a,b): get b charcters starting from index a 
                    batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                    batch_buffer.push_back(chrm_segment_seq+"\n");
                    A_segment_start = strandbreakLoc1[k]-seqStartIndex;
                    A_segment_index++; k++;
                    chrm_seg_weights.push_back(0.0);                                                            // Corresponds to the chromosome ID
                    chrm_seg_weights.push_back(static_cast<double>(A_segment_length)/total_seq_length);
                    if (GCslope != 0.0){                                                                        // Go through the calculation of GC bias only if there is a non-zero bias set
                        double chrm_GC_fraction = GCBias::get_GCfraction(chrm_segment_seq);                     // Get the GC fraction in the chromosome segment 
                        chrm_GC_bias.push_back(0.0);                                                            // Corresponds to the chromosome ID
                        chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                           // Stores the GC bias into a vector for each chromosome segment
                    }
                }else{                                                                                          // If the strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(A_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                        std::string substring = chromSeq_A.substr(A_segment_start);
                        batch_buffer.push_back(substring+"\n");
                        chrm_seg_weights.push_back(0.0);                                                        // Corresponds to the chromosome ID
                        chrm_seg_weights.push_back(static_cast<double>(substring.size())/total_seq_length);
                        if (GCslope != 0.0){                                                                    // Go through the calculation of GC bias only if there is a non-zero bias set
                            double chrm_GC_fraction = GCBias::get_GCfraction(substring);                        // Get the GC fraction in the chromosome segment 
                            chrm_GC_bias.push_back(0.0);                                                        // Corresponds to the chromosome ID
                            chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                       // Stores the GC bias into a vector for each chromosome segment
                        }
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segment tag
                        batch_buffer.push_back(chromID_A+"\n");    
                        std::string substring = chromSeq_A;
                        batch_buffer.push_back(substring+"\n");
                        chrm_seg_weights.push_back(0.0);                                                        // Corresponds to the chromosome ID
                        chrm_seg_weights.push_back(static_cast<double>(substring.size())/total_seq_length);
                        if (GCslope != 0.0){                                                                    // Go through the calculation of GC bias only if there is a non-zero bias set
                            double chrm_GC_fraction = GCBias::get_GCfraction(substring);                        // Get the GC fraction in the chromosome segment 
                            chrm_GC_bias.push_back(0.0);                                                        // Corresponds to the chromosome ID
                            chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                       // Stores the GC bias into a vector for each chromosome segment
                        }
                    }
                    break;
                }
            }    
        }else{                                                                                                  // If there are no DNA breaks in the forward strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_A+"\n"+chromSeq_A+"\n");
            chrm_seg_weights.push_back(0.0);                                                                    // Corresponds to the chromosome ID
            chrm_seg_weights.push_back(static_cast<double>(chromSeq_A.size())/total_seq_length);            
            if (GCslope != 0.0){                                                                                // Go through the calculation of GC bias only if there is a non-zero bias set
                double chrm_GC_fraction = GCBias::get_GCfraction(chromSeq_A);                                   // Get the GC fraction in the chromosome segment 
                chrm_GC_bias.push_back(0.0);                                                                    // Corresponds to the chromosome ID
                chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                                   // Stores the GC bias into a vector for each chromosome segment        
            }
        }
        
        // Processing the reverse complementary strand next 
        long B_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long B_segment_length{0};
        int B_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment in forward strand 
        if(l<strandbreakLoc2.size()){                                                                           // If there are unprocessed strand breaks present in reverse strand, then
            while(l<strandbreakLoc2.size()){                                                            
                if(strandbreakLoc2[l]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    B_segment_length = (strandbreakLoc2[l]-seqStartIndex) - B_segment_start;
                    if(B_segment_length==0){l++; continue;}                                                     // If the break is at the same location as the chromosome discontinuity, then continue to the next break
                    std::string chrm_segment_seq = chromSeq_B.substr(B_segment_start, B_segment_length);        // substr(a,b): get b charcters starting from index a 
                    batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                    batch_buffer.push_back(chrm_segment_seq+"\n");
                    B_segment_start = strandbreakLoc2[l]-seqStartIndex;                           
                    B_segment_index++; l++;
                    chrm_seg_weights.push_back(0.0);                                                            // Corresponds to the chromosome ID
                    chrm_seg_weights.push_back(static_cast<double>(B_segment_length)/total_seq_length);            
                    if (GCslope != 0.0){                                                                        // Go through the calculation of GC bias only if there is a non-zero bias set
                        double chrm_GC_fraction = GCBias::get_GCfraction(chrm_segment_seq);                     // Get the GC fraction in the chromosome segment 
                        chrm_GC_bias.push_back(0.0);                                                            // Corresponds to the chromosome ID
                        chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                           // Stores the GC bias into a vector for each chromosome segment        
                    }
                }else{                                                                                          // If strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(B_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                        std::string substring = chromSeq_B.substr(B_segment_start);
                        batch_buffer.push_back(substring+"\n");
                        chrm_seg_weights.push_back(0.0);                                                        // Corresponds to the chromosome ID
                        chrm_seg_weights.push_back(static_cast<double>(substring.size())/total_seq_length);            
                        if (GCslope != 0.0){                                                                    // Go through the calculation of GC bias only if there is a non-zero bias set
                            double chrm_GC_fraction = GCBias::get_GCfraction(substring);                        // Get the GC fraction in the chromosome segment 
                            chrm_GC_bias.push_back(0.0);                                                        // Corresponds to the chromosome ID
                            chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                       // Stores the GC bias into a vector for each chromosome segment        
                        }
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_B+"\n");
                        std::string substring = chromSeq_B;
                        batch_buffer.push_back(substring+"\n");
                        chrm_seg_weights.push_back(0.0);                                                        // Corresponds to the chromosome ID
                        chrm_seg_weights.push_back(static_cast<double>(substring.size())/total_seq_length);            
                        if (GCslope != 0.0){                                                                    // Go through the calculation of GC bias only if there is a non-zero bias set
                            double chrm_GC_fraction = GCBias::get_GCfraction(substring);                        // Get the GC fraction in the chromosome segment 
                            chrm_GC_bias.push_back(0.0);                                                        // Corresponds to the chromosome ID
                            chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                       // Stores the GC bias into a vector for each chromosome segment        
                        }
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no strand breaks in the reverse strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_B+"\n"+chromSeq_B+"\n");
            chrm_seg_weights.push_back(0.0);                                                                    // Corresponds to the chromosome ID
            chrm_seg_weights.push_back(static_cast<double>(chromSeq_B.size())/total_seq_length);            
            if (GCslope != 0.0){                                                                                // Go through the calculation of GC bias only if there is a non-zero bias set
                double chrm_GC_fraction = GCBias::get_GCfraction(chromSeq_B);                                   // Get the GC fraction in the chromosome segment 
                chrm_GC_bias.push_back(0.0);                                                                    // Corresponds to the chromosome ID
                chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                                   // Stores the GC bias into a vector for each chromosome segment        
            }
        }

        if (batch_buffer.size() >= batchSize){                                                                  // Check if the batch buffer is full, and write it to the file if needed.
            writeBatchToMMFile(batch_buffer, position_in_MM, outFileMapping, outFileSize, outFilePath);
        }
        /* // Breaking chromosomes into segments, wherever there is a DSB 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int segment_index{1};                                                                                   // Temporary index variable to name each DNA segment  
        if(DNAbreaks.size()!=0){                                                                                // If there are DNA breaks present, then
            while(i<DNAbreaks.size()){                                                            
                if(DNAbreaks[i]<=(seqStartIndex+seq_length)){                                                   // If the DNA break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                    A_segment_length = (DNAbreaks[i]-seqStartIndex) - A_segment_start;
                    outputFile<<chromSeq_A.substr(A_segment_start, A_segment_length)<<"\n";                     // substr(a,b): get b charcters starting from index a 
                    A_segment_start = DNAbreaks[i]-seqStartIndex;

                    outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                    B_segment_length = B_segment_end - (seq_length-(DNAbreaks[i]-seqStartIndex));
                    outputFile<<chromSeq_B.substr((seq_length-(DNAbreaks[i]-seqStartIndex)),B_segment_length)<<"\n";
                    B_segment_end = seq_length-(DNAbreaks[i]-seqStartIndex);
                    
                    segment_index++;
                    DNAbreaks.erase(DNAbreaks.begin() + i);
                }else{                                                                                          // If DNA break is not in the chromosome Or if it is the final segment in a chromosome
                    if(segment_index!=1){                                                                       // if it is the last segment, write ID with segment tag
                        outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }else{                                                                                      // If the DNA break is not in chromosome, write without the segmenet tag
                        outputFile<<chromID_A<<"\n";    
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no DNA breaks in the genome, then just copy the sequence as is without segmenting
            outputFile<<chromID_A<<"\n"<<chromSeq_A<<"\n"<<chromID_B<<"\n"<<chromSeq_B<<"\n";
        } */

        seqStartIndex += seq_length;
    }
    writeBatchToMMFile(batch_buffer, position_in_MM, outFileMapping, outFileSize, outFilePath);                 // If there are unwritten data in batch buffer, write that too when the loop ends
    if (msync(outFileMapping, outFileSize, MS_SYNC) == -1){                                                     // After writing your data to the memory-mapped file using mmap, before closing the file, call msync to flush/sync the changes.
        perror("\nERROR: Failed to synchronize memory-mapped data to file.\n");
    }
    munmap(outFileMapping, outFileSize);                                                                        // Unmap the memory-map to avoid memory leaks after use
    
    if(!chrm_GC_bias.empty()){                                                                                  // Modify chromosome segment weights if there is GC bias also to be considered
        double total_chrm_weight{0};                                                                            // A temporary variable to hold the total weight from each chrm weight for normalization 
        for (size_t i = 0; i<chrm_seg_weights.size(); ++i){
            chrm_seg_weights[i] *= (1+chrm_GC_bias[i]);                                                         // Length weight * GC bias weight in ragne [1,2]. GC bias is changed into this range to avoid total bias being 0 when GC bias is 0.
            total_chrm_weight += chrm_seg_weights[i];                                                           // Sum up all the individual chrm weights for normalization
        }

        for (size_t i =0; i<chrm_seg_weights.size(); i++){
            chrm_seg_weights[i] /= total_chrm_weight;                                                           // Divide each chromosome length with the total weight to normalize the vector
        }
    }
    
    return (chrm_seg_weights);
}
//--------------------------------------------------------------------------------------------
/*
    This method builds a mutated cell genome.






*/
// std::vector<double> buildMutatedCellGenome_from_MM(const std::string& outputPath, const std::string& fileName, const std::string& mutationFile, char* genomeTemplate_data, size_t templateSize, long ref_seq_length, NGSParameters& parameters){
//     size_t outFileSize = templateSize+10000;                                                                    // Temporary variable to hold the size of the current file being processed. This get dynamically changed if needed. Starts with undamagedfile size +10kb extra
//     //std::string mutFile = *parameters.get_output_directory()+"/mutations.txt";
//     std::string outFilePath = outputPath+fileName;                                                              // Path to the output fasta file that we want to write
//     char* outFileMapping = createMemoryMappedFile(outFilePath,outFileSize);                                     // Create a memory mapped output file for each damaged cell genome
//     char* position_in_MM = outFileMapping;                                                                      // Pointer to the current position in the MM as we write. Starts with the pointer to the beginning
//     std::vector<double> chrm_seg_weights;                                                                      // This is a vector to hold the weights of each chromosome segment scaled to its length
//     if (position_in_MM == nullptr){return chrm_seg_weights;}                                                    // Return if the template file memory map is not valid
//     std::vector<double> chrm_GC_bias;                                                                           // Temporary vector to hold the GC bias associated with each chromosome segment
//     long total_seq_length = ref_seq_length*2;                                                                   // Twice the length because it has two strands contributing
//     std::string chromID_A;                                                                                      // Strings to hold the chrom IDs and Sequences from the file
//     std::string chromSeq_A;
//     std::string chromID_B;
//     std::string chromSeq_B;
//     long seqStartIndex = 0;                                                                                     // Variable that will hold the starting index for each chromsome sequence in the genome
//     const int batchSize = 100;                                                                                  // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
//     std::vector<std::string> batch_buffer;                                                                      // Create a buffer for storing output data.
//     size_t position = 0;                                                                                        // Temporary variable to hold the last read position in the memory map
//     double GCslope = GCBias::get_GCbias_slope();                                                                // Temporary variable to hold the GC biase slope to avoid calling the function again and again
  
//     // std::cout << "check0" << std::endl;
//     int num_chrom = parameters.get_number_chromo();

//     std::vector<Chromosome_Metadata> cdata_matrix;

//     std::string input_mutations_path = *parameters.get_input_mutations_path();
//     std::string output_mutations_path = *parameters.get_output_directory() + "/mutations.json";

//     bool flag = check_mutation_file(input_mutations_path);

//     if (flag == true){
//         cdata_matrix = input_mutation_file(input_mutations_path).cdata_matrix;
//     }
//     else if (flag == false){
//         cdata_matrix = mutation_generator(parameters);
//         std::string cell_id = rng::random_id();
//         //std::cout << cell_id << std::endl;
//         //std::cin.get();
//         Cell_Metadata cell;
//         cell.cell_id = cell_id;
//         cell.cdata_matrix = cdata_matrix;
//         output_mutation_file(cell, output_mutations_path);
//     }
//     else{
//         //throw some error
//     }
//     std::cout << "check" << std::endl;
//     // for (int kk = 0; kk < cdata_matrix.size(); kk++){
//     //     std::cout << cdata_matrix[kk].chromosome_number << std::endl;
//     // }

//     // std::cin.get();

//     int chromCount = 1;
//     while (readFastaMemoryMap(genomeTemplate_data, templateSize, position, chromID_A, chromSeq_A, chromID_B, chromSeq_B)){// Get forward, backward sequences and their repective IDs for each chroms, one at a time
//         long seq_length = chromSeq_A.size();
//         int i = chromCount - 1;
//         cdata_matrix[i].chromA_id = chromID_A;
//         cdata_matrix[i].chromB_id = chromID_B;
//         //create mutated mirror of cell within cdata_matrix, then write it to file after the while loop is over
//         //this way I don't have to make some unreadable logic to write out to file while still mutating the cell
//         long long start_make = 0;
//         for (long long k = 0; k < seq_length; k++) { 
//             long long idx = start_make + k;
//             if (idx < 0 || idx > seq_length) {
//                 break;
//             }
//             cdata_matrix[i].chromA_seq += chromSeq_A[idx];
//         }
//         for (long long k = 0; k < seq_length; ++k) {
//             long long idx = start_make + k + 1;
//             if (idx < 0 || idx > seq_length) {
//                 break;
//             }
//             cdata_matrix[i].chromB_seq += chromSeq_B[seq_length - idx];   //front fill rev-comp leading indel chromosome
//         }
//         // for (int f = 0; f < seq_length; f++){
//         //     std::cout << cdata_matrix[i].chromA_seq[f] << cdata_matrix[i].chromB_seq[f] << " ";
//         // }
//         if (cdata_matrix[i].num_mutations == 0){//if current chromosome has no mutation
//             //do nothing lol
//         }
//         else{
//             for (int j = 0; j < cdata_matrix[i].mdata_matrix.size(); j++){
//                 //std::cout << "check4" << std::endl;
//                 if (cdata_matrix[i].mdata_matrix[j].mutation_type == "longdel"){
//                     //std::cout << "check5" << std::endl;
//                     int mutLocation = static_cast<int>(std::floor(cdata_matrix[i].mdata_matrix[j].normalized_position * static_cast<double>(seq_length)));            // pick mutation loacation on sequence
//                     cdata_matrix[i].mdata_matrix[j].position = mutLocation;
//                     int length = cdata_matrix[i].mdata_matrix[j].length;
//                     long long start = 0;
//                     cdata_matrix[i].chromA_seq.erase(mutLocation, length);
//                     cdata_matrix[i].chromB_seq.erase(mutLocation, length);
//                     seq_length = seq_length - length;
//                 }
//                 else if (cdata_matrix[i].mdata_matrix[j].mutation_type == "balinv"){
//                     //std::cout << "check6" << std::endl;
//                     std::string C1A = "";
//                     std::string C1B = "";
//                     int mutLocation = static_cast<int>(std::floor(cdata_matrix[i].mdata_matrix[j].normalized_position * static_cast<double>(seq_length)));            // pick mutation loacation on sequence
//                     cdata_matrix[i].mdata_matrix[j].position = mutLocation;
//                     int length = cdata_matrix[i].mdata_matrix[j].length;
//                     long long start_mut = mutLocation;
//                     long long start = 0;
//                     for (int k = 0; k < length; k++) {
//                         int idx = start_mut + k;
//                         if (idx < 0 || idx > seq_length) {
//                             break;
//                         }
//                         C1A += cdata_matrix[i].chromA_seq[idx];
//                     }
//                     for (int k = 0; k < length; k++) {
//                         int idx = start_mut + k;
//                         if (idx < 0 || idx > seq_length) {
//                             break;
//                         }
//                         C1B += cdata_matrix[i].chromB_seq[idx];
//                     }
//                     reverse(C1A.begin(), C1A.end());
//                     reverse(C1B.begin(), C1B.end());
//                     for (int k = 0; k < length; k++) {
//                         int idx = start_mut + k;
//                         if (idx < 0 || idx > seq_length) {
//                             break;
//                         }
//                         cdata_matrix[i].chromA_seq[idx] = C1B[k];
//                     }
//                     for (int k = 0; k < length; k++) {
//                         int idx = start_mut + k;
//                         if (idx < 0 || idx > seq_length) {
//                             break;
//                         }
//                         cdata_matrix[i].chromB_seq[idx] = C1A[k];
//                     }
//                 }
//                 else if (cdata_matrix[i].mdata_matrix[j].mutation_type == "baltrans"){
//                     //std::cout << "check7" << std::endl;
//                     std::string C1A = "";
//                     std::string C1B = "";
//                     std::string C2A = "";
//                     std::string C2B = "";
//                     if (cdata_matrix[i].chromosome_number < cdata_matrix[i].mdata_matrix[j].pair){
//                         int length = cdata_matrix[i].mdata_matrix[j].length;
//                         int mutLocation = seq_length - length - 1;            // pick mutation loacation on sequence
//                         cdata_matrix[i].mdata_matrix[j].position = mutLocation;
//                         long long start_mut = mutLocation;
//                         long long start = 0;
//                         for (long long k = mutLocation; k < seq_length; k++) {
//                             long long idx = start + k;
//                             if (idx < 0 || idx > seq_length) {
//                                 break;
//                             }
//                             C1A += cdata_matrix[i].chromA_seq[idx];
//                         }
//                         for (long long k = mutLocation + length; k < seq_length; k++){
//                             long long idx = start + k;
//                             if (idx < 0 || idx > seq_length) {
//                                 break;
//                             }
//                             C1B += cdata_matrix[i].chromB_seq[idx];
//                         }
//                         cdata_matrix[i].chromA_seq.erase(mutLocation, length);
//                         cdata_matrix[i].chromB_seq.erase(mutLocation, length);
//                         seq_length = seq_length - length + cdata_matrix[cdata_matrix[i].mdata_matrix[j].pair - 1].mdata_matrix[j].length;
//                     }
//                     else if (cdata_matrix[i].chromosome_number > cdata_matrix[i].mdata_matrix[j].pair){
//                         int length = cdata_matrix[i].mdata_matrix[j].length;
//                         int mutLocation = seq_length - length - 1;            // pick mutation loacation on sequence
//                         cdata_matrix[i].mdata_matrix[j].position = mutLocation;
//                         long long start_mut = mutLocation;
//                         long long start = 0;
//                         for (long long k = mutLocation; k < seq_length; k++) { 
//                             long long idx = start + k;
//                             if (idx < 0 || idx > seq_length) {
//                                 break;
//                             }
//                             C2A += cdata_matrix[i].chromA_seq[idx];
//                         }
//                         for (long long k = 0; k < mutLocation; k++) {
//                             long long idx = start + k;
//                             if (idx < 0 || idx > seq_length) {
//                                 break;
//                             }
//                             C2B += cdata_matrix[i].chromB_seq[idx];   //front fill rev-comp leading translocated chromosome
//                         }
//                         cdata_matrix[i].chromA_seq.erase(mutLocation, length);
//                         cdata_matrix[i].chromB_seq.erase(mutLocation, length);
//                         seq_length = seq_length - length + cdata_matrix[cdata_matrix[i].mdata_matrix[j].pair - 1].mdata_matrix[j].length;

//                         cdata_matrix[i].chromA_seq += C1A;
//                         cdata_matrix[i].chromB_seq += C1B;
//                         cdata_matrix[cdata_matrix[i].mdata_matrix[j].pair - 1].chromA_seq += C2A;
//                         cdata_matrix[cdata_matrix[i].mdata_matrix[j].pair - 1].chromA_seq += C2B;

//                     }
//                     else if (cdata_matrix[i].chromosome_number == cdata_matrix[i].mdata_matrix[j].pair){
//                         cdata_matrix[i].mdata_matrix[j].mutation_type = "none";
//                     }
//                 }
//                 else if (cdata_matrix[i].mdata_matrix[j].mutation_type == "indel"){
//                     //std::cout << "check8" << std::endl;
//                     if (cdata_matrix[i].mdata_matrix[j].inordel == "del"){
//                         int mutLocation = static_cast<int>(std::floor(cdata_matrix[i].mdata_matrix[j].normalized_position * static_cast<double>(seq_length)));            // pick mutation loacation on sequence
//                         cdata_matrix[i].mdata_matrix[j].position = mutLocation;
//                         int length = cdata_matrix[i].mdata_matrix[j].length;
//                         long long start = 0;
//                         cdata_matrix[i].chromA_seq.erase(mutLocation, length);
//                         cdata_matrix[i].chromB_seq.erase(mutLocation, length);
//                         seq_length = seq_length - length;
//                     }
//                     else if (cdata_matrix[i].mdata_matrix[j].inordel == "in"){
//                         int mutLocation = static_cast<int>(std::floor(cdata_matrix[i].mdata_matrix[j].normalized_position * static_cast<double>(seq_length)));            // pick mutation loacation on sequence
//                         cdata_matrix[i].mdata_matrix[j].position = mutLocation;
//                         int length = cdata_matrix[i].mdata_matrix[j].length;
//                         long long start = 0;
//                         const std::string random_insertion = generateDNA(length);
//                         const std::string random_insertion_forward_comp = forwardComplement(random_insertion);
//                         cdata_matrix[i].chromA_seq.insert(mutLocation, random_insertion);
//                         cdata_matrix[i].chromB_seq.insert(mutLocation, random_insertion_forward_comp);
//                     }
//                 }
//                 else if (cdata_matrix[i].mdata_matrix[j].mutation_type == "none"){
//                     continue;
//                 }
//                 else{
//                     throw std::invalid_argument("Mutation type"+cdata_matrix[i].mdata_matrix[j].mutation_type  +" unknown.");
//                 }
//             }
//         }
//         //std::cout << "check?" << std::endl;
//         chromCount++;
//         seqStartIndex += seq_length;
//     }
//     //std::cout << "check9" << std::endl;
//     for (int i = 0; i < num_chrom; i++){
//         //write-out function calls
//         //const std::string long_deletion_data = fileName + ", chromosome: " + std::to_string(chromCount) + ", long deletion at position " + std::to_string(mutLocation) + " of length " + std::to_string(length) + "\n";
//         //report_mutations(mutationFile, long_deletion_data);
//         //push mutated chromosome into buffer
//         batch_buffer.push_back(cdata_matrix[i].chromA_id+"\n"+cdata_matrix[i].chromA_seq+"\n");
//         chrm_seg_weights.push_back(0.0);                                                                    // Corresponds to the chromosome ID
//         chrm_seg_weights.push_back(static_cast<double>(cdata_matrix[i].chromA_seq.size())/total_seq_length);            
//         if (GCslope != 0.0){                                                                                // Go through the calculation of GC bias only if there is a non-zero bias set
//             double chrm_GC_fraction = GCBias::get_GCfraction(cdata_matrix[i].chromA_seq);                                   // Get the GC fraction in the chromosome segment 
//             chrm_GC_bias.push_back(0.0);                                                                    // Corresponds to the chromosome ID
//             chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                                   // Stores the GC bias into a vector for each chromosome segment        
//         }  

//         std::reverse(cdata_matrix[i].chromB_seq.begin(), cdata_matrix[i].chromB_seq.end());

//         batch_buffer.push_back(cdata_matrix[i].chromB_id+"\n"+cdata_matrix[i].chromB_seq+"\n");
//         chrm_seg_weights.push_back(0.0);                                                                    // Corresponds to the chromosome ID
//         chrm_seg_weights.push_back(static_cast<double>(cdata_matrix[i].chromB_seq.size())/total_seq_length);            
//         if (GCslope != 0.0){                                                                                // Go through the calculation of GC bias only if there is a non-zero bias set
//             double chrm_GC_fraction = GCBias::get_GCfraction(cdata_matrix[i].chromB_seq);                                   // Get the GC fraction in the chromosome segment 
//             chrm_GC_bias.push_back(0.0);                                                                    // Corresponds to the chromosome ID
//             chrm_GC_bias.push_back(GCBias::get_GCbias(chrm_GC_fraction));                                   // Stores the GC bias into a vector for each chromosome segment        
//         }
//     }
//     //std::cout << "check10" << std::endl;
//     writeBatchToMMFile(batch_buffer, position_in_MM, outFileMapping, outFileSize, outFilePath);                 // If there are unwritten data in batch buffer, write that too when the loop ends
//     if (msync(outFileMapping, outFileSize, MS_SYNC) == -1){                                                     // After writing your data to the memory-mapped file using mmap, before closing the file, call msync to flush/sync the changes.
//         perror("\nERROR: Failed to synchronize memory-mapped data to file.\n");
//     }
//     munmap(outFileMapping, outFileSize);                                                                        // Unmap the memory-map to avoid memory leaks after use
//     if(!chrm_GC_bias.empty()){                                                                                  // Modify chromosome segment weights if there is GC bias also to be considered
//         double total_chrm_weight{0};                                                                            // A temporary variable to hold the total weight from each chrm weight for normalization 
//         for (size_t i = 0; i<chrm_seg_weights.size(); ++i){
//             chrm_seg_weights[i] *= (1+chrm_GC_bias[i]);                                                         // Length weight * GC bias weight in ragne [1,2]. GC bias is changed into this range to avoid total bias being 0 when GC bias is 0.
//             total_chrm_weight += chrm_seg_weights[i];                                                           // Sum up all the individual chrm weights for normalization
//         }

//         for (size_t i =0; i<chrm_seg_weights.size(); i++){
//             chrm_seg_weights[i] /= total_chrm_weight;                                                           // Divide each chromosome length with the total weight to normalize the vector
//         }
//     }
//     return (chrm_seg_weights);
// }

namespace {

struct ReservedInterval {
    std::size_t begin = 0;
    std::size_t end = 0; // half-open interval [begin, end)
};

bool intervals_overlap(
    const ReservedInterval& first,
    const ReservedInterval& second)
{
    return first.begin < second.end && second.begin < first.end;
}

double validate_normalized_position(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Mutation position must be finite.");
    }

    return std::clamp(
        value,
        0.0,
        std::nextafter(1.0, 0.0)
    );
}

int checked_position_to_int(std::size_t position) {
    if (position > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("Mutation coordinate exceeds INT_MAX.");
    }
    return static_cast<int>(position);
}

std::size_t resolve_interval_start(
    double normalized_position,
    int requested_length,
    std::size_t sequence_length)
{
    if (requested_length <= 0) {
        throw std::invalid_argument("Mutation length must be positive.");
    }

    const std::size_t length = static_cast<std::size_t>(requested_length);
    if (length > sequence_length) {
        throw std::invalid_argument(
            "Mutation length exceeds the assigned chromosome length."
        );
    }

    const double normalized =
        validate_normalized_position(normalized_position);
    const std::size_t maximum_start = sequence_length - length;

    std::size_t start = static_cast<std::size_t>(
        std::floor(
            normalized * static_cast<double>(maximum_start + 1)
        )
    );

    return std::min(start, maximum_start);
}

std::size_t resolve_insertion_position(
    double normalized_position,
    std::size_t sequence_length)
{
    const double normalized =
        validate_normalized_position(normalized_position);

    std::size_t position = static_cast<std::size_t>(
        std::floor(
            normalized * static_cast<double>(sequence_length + 1)
        )
    );

    return std::min(position, sequence_length);
}

bool is_local_mutation(const Mutation_Metadata& mutation) {
    return mutation.mutation_type == "longdel" ||
           mutation.mutation_type == "balinv" ||
           mutation.mutation_type == "delins" ||
           mutation.mutation_type == "indel";
}

ReservedInterval proposed_interval(
    const Mutation_Metadata& mutation,
    std::size_t original_sequence_length,
    std::size_t& resolved_position)
{
    if ((mutation.mutation_type == "indel" ||
         mutation.mutation_type == "delins") &&
        mutation.inordel == "in")
    {
        resolved_position = resolve_insertion_position(
            mutation.normalized_position,
            original_sequence_length
        );

        // Reserve one anchor coordinate so an insertion is not placed inside
        // another mutation interval. The conceptual end position L is allowed.
        return {resolved_position, resolved_position + 1};
    }

    resolved_position = resolve_interval_start(
        mutation.normalized_position,
        mutation.length,
        original_sequence_length
    );

    return {
        resolved_position,
        resolved_position + static_cast<std::size_t>(mutation.length)
    };
}

void resolve_local_mutation_positions(
    Chromosome_Metadata& chromosome,
    bool allow_resampling)
{
    const std::size_t original_length = chromosome.chromA_seq.size();
    std::vector<ReservedInterval> occupied_intervals;

    for (auto& mutation : chromosome.mdata_matrix) {
        if (!is_local_mutation(mutation)) {
            continue;
        }

        constexpr int maximum_attempts = 1000;
        bool position_found = false;

        for (int attempt = 0; attempt < maximum_attempts; ++attempt) {
            std::size_t resolved_position = 0;
            const ReservedInterval candidate = proposed_interval(
                mutation,
                original_length,
                resolved_position
            );

            const bool conflict = std::any_of(
                occupied_intervals.begin(),
                occupied_intervals.end(),
                [&](const ReservedInterval& existing) {
                    return intervals_overlap(candidate, existing);
                }
            );

            if (!conflict) {
                mutation.position = checked_position_to_int(resolved_position);
                occupied_intervals.push_back(candidate);
                position_found = true;
                break;
            }

            if (!allow_resampling) {
                throw std::runtime_error(
                    "Overlapping mutations were found in the input mutation file "
                    "for chromosome " +
                    std::to_string(chromosome.chromosome_number) + "."
                );
            }

            mutation.normalized_position = rng::double_sample_0to1();
        }

        if (!position_found) {
            throw std::runtime_error(
                "Could not place a non-overlapping mutation on chromosome " +
                std::to_string(chromosome.chromosome_number) +
                " after 1000 attempts."
            );
        }
    }
}

void verify_strand_lengths(const Chromosome_Metadata& chromosome) {
    if (chromosome.chromA_seq.size() != chromosome.chromB_seq.size()) {
        throw std::runtime_error(
            "The two strands of chromosome " +
            std::to_string(chromosome.chromosome_number) +
            " have different lengths."
        );
    }
}

void apply_deletion(
    Chromosome_Metadata& chromosome,
    Mutation_Metadata& mutation)
{
    const std::size_t position = static_cast<std::size_t>(mutation.position);
    const std::size_t length = static_cast<std::size_t>(mutation.length);

    if (position > chromosome.chromA_seq.size() ||
        length > chromosome.chromA_seq.size() - position)
    {
        throw std::out_of_range("Deletion interval is outside the chromosome.");
    }

    chromosome.chromA_seq.erase(position, length);
    chromosome.chromB_seq.erase(position, length);
    verify_strand_lengths(chromosome);
}

void apply_insertion(
    Chromosome_Metadata& chromosome,
    Mutation_Metadata& mutation)
{
    const std::size_t position = static_cast<std::size_t>(mutation.position);

    if (position > chromosome.chromA_seq.size()) {
        throw std::out_of_range("Insertion position is outside the chromosome.");
    }

    std::string inserted_sequence = mutation.base_pairs;
    if (inserted_sequence.empty()) {
        inserted_sequence = generateDNA(mutation.length);
        mutation.base_pairs = inserted_sequence;
    }

    if (inserted_sequence.size() !=
        static_cast<std::size_t>(mutation.length))
    {
        throw std::runtime_error(
            "Stored insertion sequence length does not match mutation.length."
        );
    }

    const std::string complementary_sequence =
        forwardComplement(inserted_sequence);

    chromosome.chromA_seq.insert(position, inserted_sequence);
    chromosome.chromB_seq.insert(position, complementary_sequence);
    verify_strand_lengths(chromosome);
}

void apply_balanced_inversion(
    Chromosome_Metadata& chromosome,
    Mutation_Metadata& mutation)
{
    const std::size_t position = static_cast<std::size_t>(mutation.position);
    const std::size_t length = static_cast<std::size_t>(mutation.length);

    if (position > chromosome.chromA_seq.size() ||
        length > chromosome.chromA_seq.size() - position)
    {
        throw std::out_of_range("Inversion interval is outside the chromosome.");
    }

    std::string segment_a = chromosome.chromA_seq.substr(position, length);
    std::string segment_b = chromosome.chromB_seq.substr(position, length);

    std::reverse(segment_a.begin(), segment_a.end());
    std::reverse(segment_b.begin(), segment_b.end());

    chromosome.chromA_seq.replace(position, length, segment_b);
    chromosome.chromB_seq.replace(position, length, segment_a);
    verify_strand_lengths(chromosome);
}

struct DelInsRecord {
    std::size_t chromosome_index = 0;
    Mutation_Metadata* mutation = nullptr;
};

void prepare_deletion_insertions(
    std::vector<Chromosome_Metadata>& chromosomes)
{
    std::map<std::uint64_t, std::vector<DelInsRecord>> events;

    for (std::size_t chromosome_index = 0;
         chromosome_index < chromosomes.size();
         ++chromosome_index)
    {
        for (auto& mutation :
             chromosomes[chromosome_index].mdata_matrix)
        {
            if (mutation.mutation_type != "delins") {
                continue;
            }

            if (mutation.event_id == 0) {
                throw std::runtime_error(
                    "A Del-Ins record has event_id 0."
                );
            }
            if (mutation.inordel != "del" &&
                mutation.inordel != "in")
            {
                throw std::runtime_error(
                    "A Del-Ins record must have subtype 'del' or 'in'."
                );
            }

            events[mutation.event_id].push_back(
                {chromosome_index, &mutation}
            );
        }
    }

    for (auto& [event_id, records] : events) {
        if (records.size() != 2) {
            throw std::runtime_error(
                "Del-Ins event " + std::to_string(event_id) +
                " does not have exactly two records."
            );
        }

        DelInsRecord* donor = nullptr;
        DelInsRecord* recipient = nullptr;

        for (auto& record : records) {
            if (record.mutation->inordel == "del") {
                if (donor != nullptr) {
                    throw std::runtime_error(
                        "Del-Ins event " + std::to_string(event_id) +
                        " has more than one donor record."
                    );
                }
                donor = &record;
            }
            else {
                if (recipient != nullptr) {
                    throw std::runtime_error(
                        "Del-Ins event " + std::to_string(event_id) +
                        " has more than one recipient record."
                    );
                }
                recipient = &record;
            }
        }

        if (donor == nullptr || recipient == nullptr) {
            throw std::runtime_error(
                "Del-Ins event " + std::to_string(event_id) +
                " requires one donor deletion and one recipient insertion."
            );
        }

        auto& donor_chromosome =
            chromosomes.at(donor->chromosome_index);
        auto& recipient_chromosome =
            chromosomes.at(recipient->chromosome_index);
        auto& donor_mutation = *donor->mutation;
        auto& recipient_mutation = *recipient->mutation;

        if (donor_mutation.pair != recipient_chromosome.chromosome_number ||
            recipient_mutation.pair != donor_chromosome.chromosome_number)
        {
            throw std::runtime_error(
                "Del-Ins event " + std::to_string(event_id) +
                " contains non-reciprocal partner metadata."
            );
        }
        if (donor_mutation.length <= 0 ||
            donor_mutation.length != recipient_mutation.length)
        {
            throw std::runtime_error(
                "Del-Ins event " + std::to_string(event_id) +
                " has invalid or unequal donor and recipient lengths."
            );
        }

        const std::size_t donor_position =
            static_cast<std::size_t>(donor_mutation.position);
        const std::size_t length =
            static_cast<std::size_t>(donor_mutation.length);

        if (donor_position > donor_chromosome.chromA_seq.size() ||
            length > donor_chromosome.chromA_seq.size() - donor_position)
        {
            throw std::out_of_range(
                "Del-Ins donor interval is outside the chromosome."
            );
        }

        // Capture the donor segment before any chromosome is edited. This keeps
        // same-chromosome and cross-chromosome Del-Ins events in the original
        // reference-coordinate system.
        recipient_mutation.base_pairs =
            donor_chromosome.chromA_seq.substr(donor_position, length);
    }
}

void apply_local_mutations(
    std::vector<Chromosome_Metadata>& chromosomes,
    bool allow_resampling)
{
    // Resolve every local coordinate before editing any chromosome. Del-Ins
    // needs both paired records to refer to the original loaded genome.
    for (auto& chromosome : chromosomes) {
        verify_strand_lengths(chromosome);
        resolve_local_mutation_positions(chromosome, allow_resampling);
    }

    prepare_deletion_insertions(chromosomes);

    for (auto& chromosome : chromosomes) {
        std::vector<Mutation_Metadata*> local_mutations;
        for (auto& mutation : chromosome.mdata_matrix) {
            if (is_local_mutation(mutation)) {
                local_mutations.push_back(&mutation);
            }
        }

        // Coordinates were calculated against the original chromosome.
        // Applying mutations from right to left prevents coordinate shifts.
        std::sort(
            local_mutations.begin(),
            local_mutations.end(),
            [](const Mutation_Metadata* first,
               const Mutation_Metadata* second)
            {
                return first->position > second->position;
            }
        );

        for (Mutation_Metadata* mutation : local_mutations) {
            if (mutation->mutation_type == "longdel") {
                apply_deletion(chromosome, *mutation);
            }
            else if (mutation->mutation_type == "balinv") {
                apply_balanced_inversion(chromosome, *mutation);
            }
            else if ((mutation->mutation_type == "indel" ||
                      mutation->mutation_type == "delins") &&
                     mutation->inordel == "del")
            {
                apply_deletion(chromosome, *mutation);
            }
            else if ((mutation->mutation_type == "indel" ||
                      mutation->mutation_type == "delins") &&
                     mutation->inordel == "in")
            {
                apply_insertion(chromosome, *mutation);
            }
            else {
                throw std::invalid_argument(
                    "Unknown local mutation type or insertion/deletion subtype."
                );
            }
        }
    }
}

struct TranslocationRecord {
    std::size_t chromosome_index = 0;
    Mutation_Metadata* mutation = nullptr;
};

void exchange_terminal_segments(
    Chromosome_Metadata& chromosome_a,
    Mutation_Metadata& mutation_a,
    Chromosome_Metadata& chromosome_b,
    Mutation_Metadata& mutation_b)
{
    verify_strand_lengths(chromosome_a);
    verify_strand_lengths(chromosome_b);

    if (mutation_a.length <= 0 || mutation_b.length <= 0) {
        throw std::invalid_argument(
            "Translocation segment lengths must be positive."
        );
    }

    const std::size_t length_a =
        static_cast<std::size_t>(mutation_a.length);
    const std::size_t length_b =
        static_cast<std::size_t>(mutation_b.length);

    if (length_a > chromosome_a.chromA_seq.size() ||
        length_b > chromosome_b.chromA_seq.size())
    {
        throw std::invalid_argument(
            "A translocation segment is longer than its chromosome."
        );
    }

    const std::size_t total_length_before =
        chromosome_a.chromA_seq.size() + chromosome_b.chromA_seq.size();

    const std::size_t cut_a = chromosome_a.chromA_seq.size() - length_a;
    const std::size_t cut_b = chromosome_b.chromA_seq.size() - length_b;

    const std::string tail_a_forward =
        chromosome_a.chromA_seq.substr(cut_a);
    const std::string tail_a_complement =
        chromosome_a.chromB_seq.substr(cut_a);
    const std::string tail_b_forward =
        chromosome_b.chromA_seq.substr(cut_b);
    const std::string tail_b_complement =
        chromosome_b.chromB_seq.substr(cut_b);

    chromosome_a.chromA_seq.replace(cut_a, length_a, tail_b_forward);
    chromosome_a.chromB_seq.replace(cut_a, length_a, tail_b_complement);
    chromosome_b.chromA_seq.replace(cut_b, length_b, tail_a_forward);
    chromosome_b.chromB_seq.replace(cut_b, length_b, tail_a_complement);

    mutation_a.position = checked_position_to_int(cut_a);
    mutation_b.position = checked_position_to_int(cut_b);

    verify_strand_lengths(chromosome_a);
    verify_strand_lengths(chromosome_b);

    const std::size_t total_length_after =
        chromosome_a.chromA_seq.size() + chromosome_b.chromA_seq.size();

    if (total_length_before != total_length_after) {
        throw std::logic_error(
            "Balanced translocation did not conserve total sequence length."
        );
    }
}

void apply_balanced_translocations(
    std::vector<Chromosome_Metadata>& chromosomes)
{
    std::map<std::uint64_t, std::vector<TranslocationRecord>> events;

    for (std::size_t chromosome_index = 0;
         chromosome_index < chromosomes.size();
         ++chromosome_index)
    {
        for (auto& mutation :
             chromosomes[chromosome_index].mdata_matrix)
        {
            if (mutation.mutation_type != "baltrans") {
                continue;
            }

            if (mutation.event_id == 0) {
                throw std::runtime_error(
                    "A translocation record has event_id 0. Regenerate older "
                    "mutation JSON files with the revised generator."
                );
            }

            events[mutation.event_id].push_back(
                {chromosome_index, &mutation}
            );
        }
    }

    for (auto& [event_id, records] : events) {
        if (records.size() != 2) {
            throw std::runtime_error(
                "Translocation event " + std::to_string(event_id) +
                " does not have exactly two records."
            );
        }

        TranslocationRecord& first = records[0];
        TranslocationRecord& second = records[1];

        if (first.chromosome_index == second.chromosome_index) {
            throw std::runtime_error(
                "A balanced translocation cannot pair a chromosome with itself."
            );
        }

        auto& chromosome_a = chromosomes.at(first.chromosome_index);
        auto& chromosome_b = chromosomes.at(second.chromosome_index);
        auto& mutation_a = *first.mutation;
        auto& mutation_b = *second.mutation;

        if (mutation_a.pair != chromosome_b.chromosome_number ||
            mutation_b.pair != chromosome_a.chromosome_number)
        {
            throw std::runtime_error(
                "Translocation event " + std::to_string(event_id) +
                " contains non-reciprocal partner metadata."
            );
        }

        exchange_terminal_segments(
            chromosome_a,
            mutation_a,
            chromosome_b,
            mutation_b
        );
    }
}

void load_all_chromosomes(
    char* genome_template_data,
    std::size_t template_size,
    std::vector<Chromosome_Metadata>& chromosomes)
{
    std::size_t read_position = 0;
    std::size_t chromosome_index = 0;

    std::string chrom_id_a;
    std::string chrom_seq_a;
    std::string chrom_id_b;
    std::string chrom_seq_b;

    while (readFastaMemoryMap(
        genome_template_data,
        template_size,
        read_position,
        chrom_id_a,
        chrom_seq_a,
        chrom_id_b,
        chrom_seq_b))
    {
        if (chromosome_index >= chromosomes.size()) {
            throw std::runtime_error(
                "The FASTA template contains more chromosomes than "
                "number_chromo."
            );
        }

        auto& chromosome = chromosomes[chromosome_index];
        chromosome.chromA_id = chrom_id_a;
        chromosome.chromB_id = chrom_id_b;
        chromosome.chromA_seq = chrom_seq_a;

        // Store strand B in the same genomic coordinate orientation as A.
        chromosome.chromB_seq.assign(
            chrom_seq_b.rbegin(),
            chrom_seq_b.rend()
        );

        verify_strand_lengths(chromosome);
        ++chromosome_index;
    }

    if (chromosome_index != chromosomes.size()) {
        throw std::runtime_error(
            "The FASTA template contains fewer chromosomes than number_chromo."
        );
    }
}

std::string clean_fasta_id(const std::string& fasta_header) {
    std::string id = fasta_header;

    if (!id.empty() && id.front() == '>') {
        id.erase(id.begin());
    }

    const std::size_t whitespace = id.find_first_of(" \t\r\n");
    if (whitespace != std::string::npos) {
        id.resize(whitespace);
    }

    return id;
}

std::string mutation_display_name(const Mutation_Metadata& mutation) {
    if (mutation.mutation_type == "longdel") {
        return "long_deletion";
    }
    if (mutation.mutation_type == "balinv") {
        return "balanced_inversion";
    }
    if (mutation.mutation_type == "baltrans") {
        return "balanced_translocation";
    }
    if (mutation.mutation_type == "delins") {
        return "deletion_insertion";
    }
    if (mutation.mutation_type == "indel" && mutation.inordel == "in") {
        return "short_insertion";
    }
    if (mutation.mutation_type == "indel" && mutation.inordel == "del") {
        return "short_deletion";
    }

    return mutation.mutation_type;
}

std::string location_1based(const Mutation_Metadata& mutation) {
    const long long position = mutation.position;

    if ((mutation.mutation_type == "indel" ||
         mutation.mutation_type == "delins") &&
        mutation.inordel == "in")
    {
        if (position == 0) {
            return "before base 1";
        }
        return "between bases " + std::to_string(position) + " and " +
               std::to_string(position + 1);
    }

    if (mutation.mutation_type == "baltrans") {
        if (position == 0) {
            return "breakpoint before base 1";
        }
        return "breakpoint between bases " + std::to_string(position) +
               " and " + std::to_string(position + 1);
    }

    return std::to_string(position + 1) + "-" +
           std::to_string(position + mutation.length);
}

struct PartnerTruthRecord {
    int chromosome_number = 0;
    std::string contig;
    int position = -1;
};

PartnerTruthRecord find_paired_event_partner_for_report(
    const Cell_Metadata& cell,
    const Mutation_Metadata& mutation)
{
    PartnerTruthRecord result;

    const bool is_paired_event =
        mutation.mutation_type == "baltrans" ||
        mutation.mutation_type == "delins";

    if (!is_paired_event || mutation.event_id == 0) {
        return result;
    }

    for (const auto& chromosome : cell.cdata_matrix) {
        for (const auto& candidate : chromosome.mdata_matrix) {
            if (&candidate == &mutation) {
                continue;
            }
            if (candidate.mutation_type == mutation.mutation_type &&
                candidate.event_id == mutation.event_id)
            {
                result.chromosome_number = chromosome.chromosome_number;
                result.contig = clean_fasta_id(chromosome.chromA_id);
                result.position = candidate.position;
                return result;
            }
        }
    }

    return result;
}

void append_mutation_truth_report(
    const Cell_Metadata& cell,
    const std::string& mutated_fasta_name,
    const std::string& report_path)
{
    std::ofstream report(report_path, std::ios::app);
    if (!report) {
        throw std::runtime_error(
            "Could not open mutation truth report: " + report_path
        );
    }

    report << "\n# MUTATION TRUTH REPORT\n"
           << "# cell_id=" << cell.cell_id << '\n'
           << "# mutated_fasta=" << mutated_fasta_name << '\n'
           << "# Local-mutation coordinates are relative to the original "
              "forward-reference contig.\n"
           << "# Translocation coordinates are cut sites after local "
              "mutations have been applied.\n"
           << "# Del-Ins donor and recipient coordinates are relative to the "
              "original forward-reference contigs.\n"
           << "# start_0based/end_0based_exclusive use half-open intervals. "
              "Insertions have start=end.\n";

    report
        << "cell_id\tmutated_fasta\tevent_id\tchromosome_number\tcontig\t"
        << "mutation_type\tsubtype\tstart_0based\tend_0based_exclusive\t"
        << "location_1based\tlength_bp\tpartner_chromosome\tpartner_contig\t"
        << "partner_breakpoint_0based\tinserted_bases\tcoordinate_basis\n";

    std::size_t mutation_count = 0;

    std::cout << "\nExpected mutations for " << mutated_fasta_name << ":\n";

    for (const auto& chromosome : cell.cdata_matrix) {
        const std::string contig = clean_fasta_id(chromosome.chromA_id);

        for (const auto& mutation : chromosome.mdata_matrix) {
            const bool is_insertion =
                (mutation.mutation_type == "indel" ||
                 mutation.mutation_type == "delins") &&
                mutation.inordel == "in";

            const long long start_0based = mutation.position;
            const long long end_0based = is_insertion
                ? start_0based
                : start_0based + mutation.length;

            const PartnerTruthRecord partner =
                find_paired_event_partner_for_report(
                    cell,
                    mutation
                );

            const std::string event_id = mutation.event_id == 0
                ? "."
                : std::to_string(mutation.event_id);
            const std::string partner_chromosome =
                partner.chromosome_number == 0
                ? "."
                : std::to_string(partner.chromosome_number);
            const std::string partner_contig = partner.contig.empty()
                ? "."
                : partner.contig;
            const std::string partner_position = partner.position < 0
                ? "."
                : std::to_string(partner.position);
            const std::string inserted_bases = mutation.base_pairs.empty()
                ? "."
                : mutation.base_pairs;
            const std::string coordinate_basis =
                mutation.mutation_type == "baltrans"
                ? "post_local_mutation"
                : "original_reference";

            report
                << cell.cell_id << '\t'
                << mutated_fasta_name << '\t'
                << event_id << '\t'
                << chromosome.chromosome_number << '\t'
                << contig << '\t'
                << mutation_display_name(mutation) << '\t'
                << mutation.inordel << '\t'
                << start_0based << '\t'
                << end_0based << '\t'
                << location_1based(mutation) << '\t'
                << mutation.length << '\t'
                << partner_chromosome << '\t'
                << partner_contig << '\t'
                << partner_position << '\t'
                << inserted_bases << '\t'
                << coordinate_basis << '\n';

            std::cout
                << "  " << mutation_display_name(mutation)
                << " | chromosome " << chromosome.chromosome_number
                << " | " << contig
                << " | " << location_1based(mutation)
                << " | length=" << mutation.length << " bp";

            if (mutation.mutation_type == "baltrans" ||
                mutation.mutation_type == "delins")
            {
                std::cout
                    << " | partner=" << partner_contig
                    << " at 0-based breakpoint " << partner_position;
            }

            std::cout << '\n';
            ++mutation_count;
        }
    }

    report << "# mutation_records=" << mutation_count << "\n";

    std::cout << "Mutation truth appended to: " << report_path << "\n";
}

void write_metadata_without_genome_sequences(
    const Cell_Metadata& cell,
    const std::string& output_path)
{
    Cell_Metadata metadata_only = cell;

    for (auto& chromosome : metadata_only.cdata_matrix) {
        chromosome.chromA_seq.clear();
        chromosome.chromB_seq.clear();
        chromosome.num_mutations =
            static_cast<int>(chromosome.mdata_matrix.size());
    }

    // The existing API takes a non-const reference.
    output_mutation_file(metadata_only, output_path);
}

std::vector<double> write_mutated_fasta_and_get_weights(
    const std::string& output_file_path,
    const std::vector<Chromosome_Metadata>& chromosomes)
{
    std::ofstream output_file(
        output_file_path,
        std::ios::binary | std::ios::trunc
    );

    if (!output_file) {
        throw std::runtime_error(
            "Could not open mutated FASTA output file: " + output_file_path
        );
    }

    std::vector<double> weights;
    std::vector<double> gc_bias;
    weights.reserve(chromosomes.size() * 4);
    gc_bias.reserve(chromosomes.size() * 4);

    const bool use_gc_bias = GCBias::get_GCbias_slope() != 0.0;

    for (const auto& chromosome : chromosomes) {
        verify_strand_lengths(chromosome);

        // Convert aligned strand B back to its FASTA 5'-to-3' orientation.
        const std::string output_strand_b(
            chromosome.chromB_seq.rbegin(),
            chromosome.chromB_seq.rend()
        );

        output_file << chromosome.chromA_id << '\n'
                    << chromosome.chromA_seq << '\n'
                    << chromosome.chromB_id << '\n'
                    << output_strand_b << '\n';

        if (!output_file) {
            throw std::runtime_error(
                "Failed while writing mutated FASTA file: " + output_file_path
            );
        }

        weights.push_back(0.0);
        weights.push_back(
            static_cast<double>(chromosome.chromA_seq.size())
        );
        weights.push_back(0.0);
        weights.push_back(
            static_cast<double>(output_strand_b.size())
        );

        if (use_gc_bias) {
            gc_bias.push_back(0.0);
            gc_bias.push_back(
                GCBias::get_GCbias(
                    GCBias::get_GCfraction(chromosome.chromA_seq)
                )
            );
            gc_bias.push_back(0.0);
            gc_bias.push_back(
                GCBias::get_GCbias(
                    GCBias::get_GCfraction(output_strand_b)
                )
            );
        }
    }

    if (use_gc_bias) {
        for (std::size_t i = 0; i < weights.size(); ++i) {
            weights[i] *= 1.0 + gc_bias[i];
        }
    }

    const double total_weight =
        std::accumulate(weights.begin(), weights.end(), 0.0);

    if (!(total_weight > 0.0)) {
        throw std::runtime_error("Mutated genome has zero total sequence length.");
    }

    for (double& weight : weights) {
        weight /= total_weight;
    }

    return weights;
}

} // namespace

std::vector<double> buildMutatedCellGenome_from_MM(
    const std::string& outputPath,
    const std::string& fileName,
    const std::string& mutationFile,
    char* genomeTemplate_data,
    std::size_t templateSize,
    long ref_seq_length,
    NGSParameters& parameters)
{
    // mutationFile is the human-readable truth-report path. Mutation JSON
    // input still comes from input_mutations_path in the parameter file.
    (void)ref_seq_length;

    if (genomeTemplate_data == nullptr) {
        throw std::invalid_argument("genomeTemplate_data is null.");
    }

    const int number_of_chromosomes = parameters.get_number_chromo();
    if (number_of_chromosomes <= 0) {
        throw std::invalid_argument("number_chromo must be positive.");
    }

    const std::string input_mutations_path =
        *parameters.get_input_mutations_path();
    const std::string cell_name =
        std::filesystem::path(fileName).stem().string();
    const std::string output_mutations_path =
        *parameters.get_output_directory() + "/" +
        cell_name + "_mutations.json";
    const std::string output_fasta_path = outputPath + fileName;

    const bool loaded_existing_metadata =
        check_mutation_file(input_mutations_path);

    Cell_Metadata cell;

    if (loaded_existing_metadata) {
        cell = input_mutation_file(input_mutations_path);
    }
    else {
        cell.cell_id = rng::random_id();
        cell.cdata_matrix = mutation_generator(parameters);
    }

    if (cell.cdata_matrix.size() !=
        static_cast<std::size_t>(number_of_chromosomes))
    {
        throw std::runtime_error(
            "Mutation metadata chromosome count does not match number_chromo."
        );
    }

    load_all_chromosomes(
        genomeTemplate_data,
        templateSize,
        cell.cdata_matrix
    );

    // Generated mutations may be resampled to avoid overlap. Positions from an
    // explicitly supplied JSON file are treated as fixed and overlap is an error.
    apply_local_mutations(
        cell.cdata_matrix,
        !loaded_existing_metadata
    );

    // Balanced translocations are applied after local mutations, including
    // Del-Ins donor deletions and recipient insertions.
    apply_balanced_translocations(cell.cdata_matrix);

    // At this point every mutation has a final absolute position, and
    // insertion sequences and translocation partner breakpoints are known.
    append_mutation_truth_report(
        cell,
        fileName,
        mutationFile
    );

    std::vector<double> weights =
        write_mutated_fasta_and_get_weights(
            output_fasta_path,
            cell.cdata_matrix
        );

    // Save final absolute positions and generated insertion bases, but omit the
    // full chromosome sequences from the JSON file.
    write_metadata_without_genome_sequences(
        cell,
        output_mutations_path
    );

    return weights;
}


std::string generateDNA(int n) { //generate random string for insertions
    const char bases[] = {'A', 'T', 'C', 'G'};
    std::string dna;

    static bool seeded = false;
    if (!seeded) {
        std::srand(std::time(nullptr));
        seeded = true;
    }

    for (int i = 0; i < n; ++i) {
        dna += bases[std::rand() % 4];
    }

    return dna;
}

std::string reverseComplement(const std::string& dna) { //generate rev-comp of a string
    std::string result;
    result.reserve(dna.size()); // avoid reallocations

    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        switch (*it) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default:  result += 'N'; 
        }
    }

    return result;
}

std::string forwardComplement(const std::string& dna) { //generate rev-comp of a string
    std::string result;
    result.reserve(dna.size()); // avoid reallocations

    for (char c : dna) {
        switch (c) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default:  result += 'N';
        }
    }

    return result;
}