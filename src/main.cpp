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
#include <iostream>
#include <fstream>
#include <string>
#include "HaplotypeVcfParser.h"



int BLOCK_SIZE { 64 };

int main(int argc, char* argv[])
{
    if (argc != 2)
        throw("Must specify vcf");

    // open VCF file and parse meta data and header
    HaplotypeVcfParser vcf_data { argv[1] };


    int n_blocks { 1 + static_cast<int>(vcf_data.n_samples() / BLOCK_SIZE) };


    // instantiate matrices to hold calculations
    Matrix covariance { vcf_data.n_samples(), vcf_data.n_samples() };

    // instantiate record object
    HaplotypeDataRecord record { vcf_data.n_samples(), vcf_data.k_founders() };

    int ii,jj,i,j,end_block_i,end_block_j,start_block_i,start_block_j;
    while(vcf_data.load_record(record)) {

        // for each founder, compute first and second moments

        for (i = 0; i < n_blocks; i++) {
            for (j = i; j < n_blocks; j++) {

                start_block_i = i*BLOCK_SIZE;
                end_block_i = (i + 1)*BLOCK_SIZE;
                if (end_block_i > vcf_data.n_samples())
                    end_block_i = vcf_data.n_samples();
                
                start_block_j = j*BLOCK_SIZE;
                end_block_j = (i + 1)*BLOCK_SIZE;
                if (end_block_j > vcf_data.n_samples())
                    end_block_j = vcf_data.n_samples();

                for (ii = start_block_i; ii < end_block_i; ii++)
                    for (jj=start_block_j; jj < end_block_j; jj++) {

                        for (int k = 0; k < vcf_data.k_founders(); k++)
                            covariance(ii, jj) += record(ii, k) * record(jj, k);
                    }

                if (i != j)
                    covariance(jj, ii) = covariance(ii,jj);
            }
        }

    }


    for (int i = 0; i < 10; i++) {
        std::cout << std::endl;

        for (int j = 0; j < 10; j++)
            std::cout << covariance(i, j) << ", ";

    }
    return 0;
}
