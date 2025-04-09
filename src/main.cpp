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
#include <cstdio>
#include <chrono>
#include "HaplotypeVcfParser.h"



int BLOCK_SIZE { 64 };
size_t MARKER_PRINT_INTERVAL { 0 };
char HELP_LONG_FLAG[] { "--help" };
char HELP_SHORT_FLAG[] { "-h" };

int main(int argc, char* argv[])
{

    if (argc != 2 && argc != 3)
        throw("Must specify vcf");

    if (argc == 2 
            && (strcmp(argv[1], HELP_SHORT_FLAG) == 0
                || strcmp(argv[1], HELP_LONG_FLAG) == 0)) {
        printf("hgrm - Compute GRM from expected haplotype counts.\n"
               "Usage\n"
               "\n"
               "  hgrm <input_vcf_filename> [<output_vcf_filename>]\n"
               "\n"
               "Options\n"
               "  output_vcf_filename   Filename to print covariance matrix\n"
               "\n"
               "Description\n"
               "  A program to compute a genetic relationship matrix from a vcf\n"
               "  with expected haplotype counts record per sample per locus.\n");


        return 0;
    }

    char* filename_input { argv[1] };
    char* filename_output { nullptr };

    if(argc == 3)
        filename_output = argv[2];


    const std::chrono::time_point timer
    { std::chrono::steady_clock::now() };


    fprintf(stdout, "Allocating memory\n");

    // open VCF file and parse meta data and header
    HaplotypeVcfParser vcf_data { filename_input, 100000 };


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

    fprintf(stdout, "Computing matrix\n");

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

        if (m_markers % MARKER_PRINT_INTERVAL == 0)
            fprintf(stdout, "Completed %zu marker loci\n", m_markers);

        m_markers++;

    }


    FILE* fout = stdout;

    if (argc == 3 && filename_output != nullptr) {

        fout = fopen(filename_output, "w");
        fprintf(stdout, "Writing results to file %s\n", filename_output);

    } else if (argc == 3 && filename_output == nullptr)
        throw std::runtime_error("Output filename is not specified");


    int i { 0 };
    int j { 0 };
    for (i = 0; i < n_samples; i++) {

        for (j = 0; j < n_samples-1; j++)
            fprintf(fout, "%0.5f,", covariance(i,j));

        fprintf(fout,"%0.5f\n", covariance(i, j));
    }


    std::chrono::steady_clock::duration delta_t
        { std::chrono::steady_clock::now() - timer };

    fprintf(stdout, "Elapsed time: %lld seconds\n",
            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());

    return 0;
}
