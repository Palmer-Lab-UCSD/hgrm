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

    Matrix first_moment { vcf_data.n_samples(), vcf_data.k_founders() };
    Matrix delta { vcf_data.n_samples(), vcf_data.k_founders() };


    // instantiate record object
    HaplotypeDataRecord record { vcf_data.n_samples(), vcf_data.k_founders() };

    // analyze each line, i.e. position, in the VCF
    size_t m_markers { 1 };

    while(vcf_data.load_record(record)) {

        // for each founder, compute first and second moments
        
        if (m_markers == 1) {
            for (int i = 0; i < vcf_data.n_samples(); i++) 
                for (int k = 0; k < vcf_data.k_founders(); k++) {
                    first_moment(i, k) = record(i, k);
                    delta(i, k) = 0;
                } 
            m_markers++;
            continue;
        } else {

            for (int i = 0; i < vcf_data.n_samples(); i++) {
                for (int k = 0; k < vcf_data.k_founders(); k++) {
                    delta(i, k) = record(i, k) - first_moment(i, k);
                    first_moment(i, k) += delta(i, k) / m_markers;
                }
            }
        }


        for (int i = 0; i < vcf_data.n_samples(); i++) {
            for (int j = i; j < vcf_data.n_samples(); j++) {
                for (int k = 0; k < vcf_data.k_founders(); k++)
                    covariance(i, j) = (m_markers-2)*covariance(i, j) / (m_markers-1) 
                                        + delta(i,k)*delta(j,k)/m_markers;

                if (i != j)
                    covariance(j, i) = covariance(i,j);
            }
        }

        m_markers++;
    }


    for (int i = 0; i < vcf_data.n_samples(); i++) {

        std::cout << std::endl;

        for (int j = 0; j < vcf_data.n_samples(); j++)
            std::cout << covariance(i, j) << ", ";

    }
    return 0;
}
