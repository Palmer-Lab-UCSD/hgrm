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




int main(int argc, char* argv[])
{
    if (argc != 2)
        throw("Must specify vcf");

    // open VCF file and parse meta data and header
    HaplotypeVcfParser vcf_data { argv[1] };

    // instantiate matrices to hold calculations

    Matrix covariance { vcf_data.n_samples(), vcf_data.n_samples() };

    Matrix first_moment { vcf_data.k_founders(), vcf_data.n_samples() };
    Matrix delta { vcf_data.k_founders(), vcf_data.n_samples() };


    // instantiate record object
    HaplotypeDataRecord record { vcf_data.k_founders(), vcf_data.n_samples() };

    // analyze each line, i.e. position, in the VCF
    size_t m_markers { 1 };

    while(vcf_data.load_record(record)) {

        // for each founder, compute first and second moments
        for (int k = 0; k < vcf_data.k_founders(); k++) {

            
            if (m_markers == 1) {
                for (int i = 0; i < vcf_data.n_samples(); i++) {
                    first_moment(k, i) = record(k, i);
                    delta(k, i) = 0;
                } 
                continue;
            }


            for (int i = 0; i < vcf_data.n_samples(); i++) {

                delta(k,i) = record(k, i) - first_moment(k, i);

                for (int j = i; j < vcf_data.n_samples(); j++) {

                    delta(k,j) = record(k, j) - first_moment(k, j);

                    covariance(i, j) = (m_markers-2) * covariance(i, j) / (m_markers-1) 
                                            + delta(k,i)*delta(k,j)/m_markers;

                    if (i != j)
                        covariance(j, i) = covariance(i,j);
                }

                first_moment(k, i) += delta(k, i) / m_markers;
            }
        }

        m_markers++;
    }


    for (int i = 0; i < vcf_data.n_samples(); i++) {

        std::cout << std::endl;

        for (int j = 0; j < vcf_data.n_samples(); j++)
            std::cout << covariance(i, j) << ", " << std::endl;

    }
    return 0;
}
