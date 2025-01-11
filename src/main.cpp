// Compute the genomic relationship matrix using haplotypes
//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
// Input argument
//    filename: vcf with haplotpye
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
    std::array<Matrix, vcf_data.k_founders()> second_moment;

    for(int k = 0; k < vcf_data.k_founders(); k++)
        second_moment[k] { vcf_data.n_samples(), vcf_data.n_samples() };



    // instantiate record object
    HaplotypeDataRecord record {};

    // analyze each line, i.e. position, in the VCF
    size_t m_markers { 0 };

    while(vcf_data.get_record(record)) {

        // for each founder, compute first and second moments
        for (int k = 0; k < vcf_data.k_founders(); k++)

            for (int i = 0; i < vcf_data.n_samples(); i++) {

                first_moment(k, i) += record(k, i);
                second_moment[k](i, i) += record(k, i) * record(k,i);

                for (int j = i+1; j < vcf_data.k_founders(); j++) {
                    second_moment[k](i, j) += record(k, i) * record(k, j);
                }
            }

        m_markers++;
    }

    // Use the first and second moments to compute the unbiased covariance 

    for (int k = 0; k < vcf_data.k_founders(); k++) {

        for (int i=0; i < vcf_data.n_samples(); i++) {
            first_moment(k, i) =  first_moment(k,i) / m_markers;

            covariance[i, i] += (second_moment[k](i,i)
                                    - first_moment(k, i)^2) / (m_markers-1);

            std::cout << covariance[i,i] << std::endl;

            for(int j=i+1; j < vcf_data.n_samples(); j++) {
                covariance[i, j] += (second_moment[k](i, j)
                                    - first_moment(k,i)*first_moment(k,j)) / (m_markers-1);

                covariance[j, i] = covariance[i, j];

                std::cout << ", " << covariance[j,i];
            }
            std::cout << std::endl;
        }
    }


    return 0;
}
