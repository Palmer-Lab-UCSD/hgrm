// Compute the genomic relationship matrix using haplotypes
//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
// Input argument
//    filename: vcf with haplotpye
//
// This program performs a single-pass computation of the 
// haplotype based genomic relationship matrix.  The approach
// is well defined for the covariance, however under my definition
// of the haplotype based covariance I had to derive the recursion
// relations myself.
//
//
//
// Acknowledgment
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
//
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include "HaplotypeVcfParser.h"




int main(int argc, char* argv[])
{
    if (argc != 2)
        throw("Must specify vcf");

    // open VCF file and parse meta data and header
    HaplotypeVcfParser vcf_data { argv[1] };

    // instantiate matrices to hold calculations

    Matrix covariance { vcf_data.n_samples(), vcf_data.n_samples() };


    #pragma omp parallel
    {
        // instantiate record object
        HaplotypeDataRecord record { vcf_data.n_samples(), vcf_data.k_founders() };

        // std::cout << std::thread::hardware_concurrency()<< std::endl;
        
    // mutex on reading records, guarantees that each thread has a unique line
    #pragma omp critical
    bool record_read { vcf_data.load_record(record) };

    while(record_read) {

        // for each founder, compute first and second moments
        //auto start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < vcf_data.n_samples(); i++) {
            for (int j = i; j < vcf_data.n_samples(); j++) {
                for (int k = 0; k < vcf_data.k_founders(); k++)
                    covariance(i, j) += record(i, k) * record(j, k);

                if (i != j)
                    covariance(j, i) = covariance(i,j);
            }
        }

        //auto end = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
       // std::cout << "Computation time: " << duration.count() << "s\n";

    }

    }

    for (int i = 0; i < vcf_data.n_samples(); i++) {

        std::cout << std::endl;

        for (int j = 0; j < vcf_data.n_samples(); j++)
            std::cout << covariance(i, j) << ", ";

    }
    return 0;
}
