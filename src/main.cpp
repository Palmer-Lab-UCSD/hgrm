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
#include <argparse.h>
#include <cstdio>
// #include <chrono>
#include <cstdlib>
#include <optional>
#include <string>


size_t MARKER_PRINT_INTERVAL { 1000 };
char HELP_LONG_FLAG[] { "--help" };
char HELP_SHORT_FLAG[] { "-h" };

int main(int argc, char* argv[])
{
    if (argc != 2 && argc != 4) {
        fprintf(stderr, "Incorrect input, see --help for correct usage.\n");
        exit(EXIT_FAILURE);
    }

    argparse::ArgParser parser {
        "hgrm: Haplotype Genetic Relationship Matrix",
        "This program computes the haplotype genetic relationship matrix"
        "from the expected haplotype counts per locus per sample and stored"
        "as a text file in the variant call format (VCF)."
    };

    parser.add_arg("--sample_names",
            argparse::ArgType::STRING,
            "The path and name of the file containing sample names to be"
            " included in computing the relationship matrix.  The file must"
            " include a single sample filename, and if necessary file system"
            " path, per line.");
    parser.add_arg("-o", 
            argparse::ArgType::STRING,
            "the path and filename that the resulting haplotype genetic"
            "relationship matrix is printed.");
    parser.add_arg("vcf",
            argparse::ArgType::STRING, 
            "the path and filename of the vcf in which the hgrm is computed.");

    if (parser.parse_args(argc, argv) != argparse::ArgStatus::SUCCESS) {
        fprintf(stderr, "Error: couldn't parse command line args, exiting\n");
        exit(EXIT_FAILURE);
    }

    std::optional<std::string> tmp {};
    if((tmp = parser.get<std::string>("vcf")) == std::nullopt) {
        fprintf(stderr, "Error retrieving vcf name");
        exit(EXIT_FAILURE);
    }
    std::string vcf_fname { tmp.value() };

    if ((tmp = parser.get<std::string>("o")) == std::nullopt) {
        fprintf(stderr, "Error retrieving output name");
        exit(EXIT_FAILURE);
    }
    std::string out_fname { tmp.value() };

    if (out_fname.size() == 0)
        out_fname = vcf_fname + ".mat";

    std::string samp_fname {};
    if ((tmp = parser.get<std::string>("sample_names")) == std::nullopt) {
        fprintf(stderr, "Error retrieving sample_names file.\n");
        exit(EXIT_FAILURE);
    }
    samp_fname = tmp.value();


    fprintf(stdout, "BCF/VCF file name: %s\n", vcf_fname.c_str());
    if (samp_fname.size() == 0)
        fprintf(stdout, "Sample file: None, use all samples\n");
    else
        fprintf(stdout, "Sample file: %s\n", samp_fname.c_str());
    fprintf(stdout, "Output matrix file: %s\n", out_fname.c_str());

    // const std::chrono::time_point timer;
    // { std::chrono::steady_clock::now() };
    
//    HaplotypeVcfParser vcf_data { filename_input, 100000 };

    fprintf(stdout, "Allocating memory\n");
    // instantiate matrices to hold calculations
//     Matrix covariance { vcf_data.n_samples(), vcf_data.n_samples() };

//
//    // open VCF file and parse meta data and header
//    HaplotypeVcfParser vcf_data { filename_input, 100000 };
//
//

//    // instantiate record object
//    HaplotypeDataRecord record { vcf_data.n_samples(), vcf_data.k_founders() };
//
//    // analyze each line, i.e. position, in the VCF
//    size_t m_markers { 1 };
//
//    double sum { 0 };
//    const double* rowi { nullptr };
//    const double* rowj { nullptr };
//    double* rowi_cov { nullptr };
//    const size_t k_founders { vcf_data.k_founders() };
//    const size_t n_samples { vcf_data.n_samples() };
//
//    std::chrono::steady_clock::duration delta_t
//        { std::chrono::steady_clock::now() - timer };
//
//    fprintf(stdout, "Computing matrix, elapsed time %lld second(s)\n",
//            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//
//    while(vcf_data.load_record(record)) {
//
//        // for each founder, compute first and second moments
//        for (size_t i = 0; i < n_samples; i++) {
//
//            rowi = &record(i, 0);
//            rowi_cov = &covariance(i, 0);
//
//            for (size_t j = i; j < n_samples; j++) {
//
//                rowj = &record(j,0);
//                sum = 0;
//
//                for (int k = 0; k < k_founders; k++)
//                    sum += rowi[k] * rowj[k];
//
//                rowi_cov[j] += sum;
//            }
//        }
//
//        if (m_markers % MARKER_PRINT_INTERVAL == 0) {
//            delta_t = std::chrono::steady_clock::now() - timer;
//
//            fprintf(stdout, "Completed %zu marker loci, elapsed time %lld second(s)\n",
//                    m_markers,
//                    std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//        }
//
//        m_markers++;
//
//    }
//
//
//    FILE* fout = stdout;
//
//    if (argc == 3 && filename_output != nullptr) {
//
//        if ((fout = fopen(filename_output, "w")) == nullptr)
//            throw std::runtime_error("Error in opening file for writing.");
//
//        delta_t = std::chrono::steady_clock::now() - timer;
//        fprintf(stdout, "Writing results to file %s, elapsed time %lld second(s)\n",
//                filename_output,
//                std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//
//    } else if (argc == 3 && filename_output == nullptr)
//        throw std::runtime_error("Output filename is not specified");
//
//
//    size_t i { 0 };
//    size_t j { 0 };
//    for (i = 0; i < n_samples; i++) {
//
//        for (j = 0; j < n_samples-1; j++) {
//            if (j < i)
//                fprintf(fout, "%0.5f,", covariance(j,i));
//            else
//                fprintf(fout, "%0.5f,", covariance(i,j));
//
//        }
//
//        fprintf(fout,"%0.5f\n", covariance(i, j));
//    }
//
//    fclose(fout);
//
//
//    delta_t = std::chrono::steady_clock::now() - timer;
//
//    fprintf(stdout, "Done, elapsed time %lld second(s)\n",
//            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//
    return 0;
}
