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
size_t MARKER_PRINT_INTERVAL { 1000 };
char HELP_LONG_FLAG[] { "--help" };
char HELP_SHORT_FLAG[] { "-h" };

int main(int argc, char* argv[])
{

    if (argc != 2 && argc != 3)
        throw std::runtime_error("Must specify vcf");

    if (argc == 2 
            && (strcmp(argv[1], HELP_SHORT_FLAG) == 0
                || strcmp(argv[1], HELP_LONG_FLAG) == 0)) {
        printf("hgrm - Compute GRM from expected haplotype counts.\n"
               "Usage\n"
               "\n"
               "  hgrm <input_vcf_filename> [<output_matrix_filename>]\n"
               "\n"
               "Options\n"
               "  output_matrix_filename   Filename to print covariance matrix\n"
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


    // instantiate matrices to hold calculations
    Matrix covariance { vcf_data.n_samples(), vcf_data.n_samples() };

    // instantiate record object
    HaplotypeDataRecord record { vcf_data.n_samples(), vcf_data.k_founders() };

    // analyze each line, i.e. position, in the VCF
    size_t m_markers { 1 };

    double sum { 0 };
    const double* rowi { nullptr };
    const double* rowj { nullptr };
    double* rowi_cov { nullptr };
    const size_t k_founders { vcf_data.k_founders() };
    const size_t n_samples { vcf_data.n_samples() };

    std::chrono::steady_clock::duration delta_t
        { std::chrono::steady_clock::now() - timer };

    fprintf(stdout, "Computing matrix, elapsed time %lld second(s)\n",
            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());

    while(vcf_data.load_record(record)) {

        // for each founder, compute first and second moments
        for (int i = 0; i < n_samples; i++) {

            rowi = &record(i, 0);
            rowi_cov = &covariance(i, 0);

            for (int j = i; j < n_samples; j++) {

                rowj = &record(j,0);
                sum = 0;

                for (int k = 0; k < k_founders; k++)
                    sum += rowi[k] * rowj[k];

                rowi_cov[j] += sum;
            }
        }

        if (m_markers % MARKER_PRINT_INTERVAL == 0) {
            delta_t = std::chrono::steady_clock::now() - timer;

            fprintf(stdout, "Completed %zu marker loci, elapsed time %lld second(s)\n",
                    m_markers,
                    std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
        }

        m_markers++;

    }


    FILE* fout = stdout;

    if (argc == 3 && filename_output != nullptr) {

        if ((fout = fopen(filename_output, "w")) == nullptr)
            throw std::runtime_error("Error in opening file for writing.");

        delta_t = std::chrono::steady_clock::now() - timer;
        fprintf(stdout, "Writing results to file %s, elapsed time %lld second(s)\n",
                filename_output,
                std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());

    } else if (argc == 3 && filename_output == nullptr)
        throw std::runtime_error("Output filename is not specified");


    size_t i { 0 };
    size_t j { 0 };
    for (i = 0; i < n_samples; i++) {

        for (j = 0; j < n_samples-1; j++) {
            if (j < i)
                fprintf(fout, "%0.5f,", covariance(j,i));
            else
                fprintf(fout, "%0.5f,", covariance(i,j));

        }

        fprintf(fout,"%0.5f\n", covariance(i, j));
    }

    fclose(fout);


    delta_t = std::chrono::steady_clock::now() - timer;

    fprintf(stdout, "Done, elapsed time %lld second(s)\n",
            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());

    return 0;
}
