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
    HaplotypeVcfParser vcf_data { argv[1], 100000 };


    int n_blocks { 1 + static_cast<int>(vcf_data.n_samples() / BLOCK_SIZE) };


    // instantiate matrices to hold calculations
    Matrix covariance { vcf_data.n_samples(), vcf_data.n_samples() };

    // instantiate record object
    HaplotypeDataRecord record { vcf_data.n_samples(), vcf_data.k_founders() };

    // analyze each line, i.e. position, in the VCF
    size_t m_markers { 1 };

    double sum { 0 };
    const double* rowi { nullptr };
    const double* rowj { nullptr };
    const size_t k_founders { vcf_data.k_founders() };
    const size_t n_samples { vcf_data.n_samples() };

    while(vcf_data.load_record(record)) {

        // for each founder, compute first and second moments
        for (int i = 0; i < n_samples; i++) {

            rowi = &record(i, 0);

            for (int j = i; j < n_samples; j++) {

                rowj = &record(j,0);
                sum = 0;

                for (int k = 0; k < k_founders; k++)
                    sum += rowi[k] * rowj[k];

                covariance(i, j) += sum;

                if (i != j)
                    covariance(j, i) += sum;
            }
        }

    }


    int i { 0 };
    int j { 0 };
    for (i = 0; i < n_samples; i++) {

        for (j = 0; j < n_samples-1; j++)
            std::cout << covariance(i, j) << ",";

        std::cout << covariance(i,j) << std::endl;
    }
    return 0;
}
