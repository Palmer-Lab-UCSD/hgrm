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
#include <argparse.h>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>

#include <logger.h>
#include <calc.h>
#include <matrix.h>


#define FAILED_CALC -1
#define SUCCESS_CALC 0


const size_t STR_BUF_LEN { 500 };
char STR_BUF[STR_BUF_LEN];


int main(int argc, char* argv[])
{
    // if (argc != 2 && argc != 4) {
    //     fprintf(stderr, "Incorrect input, see --help for correct usage.\n");
    //     exit(EXIT_FAILURE);
    // }

    argparse::ArgParser parser {
        "hgrm: Haplotype Genetic Relationship Matrix",
        "This program computes the haplotype genetic relationship matrix"
        " from the expected haplotype counts per locus per sample and stored"
        " as a text file in the variant call format (VCF)."
    };

    parser.add_arg("-o", 
            argparse::ArgType::STRING,
            "the path and filename that the resulting haplotype genetic"
            " relationship matrix is printed.");

    parser.add_arg("--sample_names",
            argparse::ArgType::STRING,
            "The path and name of the file containing sample names to be"
            " included in computing the relationship matrix.  The file must"
            " include a single sample filename, and if necessary file system"
            " path, per line.");

    parser.add_arg("--gt",
            argparse::ArgType::BOOLEAN,
            "Use sample genotypes to compute the relationship matrix");

    parser.add_arg("--eac",
            argparse::ArgType::BOOLEAN,
            "Use sample expected alt allele count to compute the"
            " relationship matrix");

    parser.add_arg("-b",
            argparse::ArgType::BOOLEAN,
            "Use both the expected alternative allele and haplotype counts to"
            " compute relationship matrix");

    parser.add_arg("--loco",
            argparse::ArgType::STRING,
            "Directory with chromosome matrix files to compute the"
            " leave-one-chromosome-out (LOCO) relationship matrix.")

    parser.add_arg("--vcf",
            argparse::ArgType::STRING, 
            "the path and filename of the vcf in which the hgrm is computed.");

    if (parser.parse_args(argc, argv) != argparse::ArgStatus::SUCCESS) {
        fprintf(stderr, "Error: couldn't parse command line args, exiting\n");
        exit(EXIT_FAILURE);
    }

    // TODO: Update below to use logger
    //
    std::optional<std::string> tmp_str {};
    if((tmp_str = parser.get<std::string>("vcf")) == std::nullopt) {
        fprintf(stderr, "Error retrieving vcf name");
        exit(EXIT_FAILURE);
    }
    std::string vcf_fname { tmp_str.value() };

    if ((tmp_str = parser.get<std::string>("o")) == std::nullopt) {
        fprintf(stderr, "Error retrieving output name");
        exit(EXIT_FAILURE);
    }
    std::string out_fname { tmp_str.value() };

    if (out_fname.size() == 0)
        out_fname = vcf_fname + ".mat";

    if ((tmp_str = parser.get<std::string>("sample_names")) == std::nullopt) {
        fprintf(stderr, "Error retrieving sample_names file.\n");
        exit(EXIT_FAILURE);
    }
    std::string samp_fname { tmp_str.value() };

    
    std::optional<bool> tmp_bool {};
    if ((tmp_bool = parser.get<bool>("gt")) == std::nullopt) {
        fprintf(stderr, "Error retrieving relationship matrix type.\n");
        exit(EXIT_FAILURE);
    }
    bool use_gt { tmp_bool.value() };

    if ((tmp_bool = parser.get<bool>("b")) == std::nullopt) {
        fprintf(stderr, "Error retrieving relationship matrix type.\n");
        exit(EXIT_FAILURE);
    }
    bool use_both { tmp_bool.value() };

    if ((tmp_bool = parser.get<bool>("eac")) == std::nullopt) {
        fprintf(stderr, "Error retrieving relationship matrix type.\n");
        exit(EXIT_FAILURE);
    }
    bool use_eac { tmp_bool.value() };


    if ((use_gt && use_both) || (use_gt && use_ds) || (use_both && use_ds)) {
        fprintf(stderr, "user must specify either use_gt, use_both, use_ds,"
                " or omit both options to compute the haplotype based"
                " relationship matrix.");
        exit(EXIT_FAILURE);
    }


    Logger log {};
    
    log.info("BCF/VCF file name: %s", vcf_fname.c_str());
    if (samp_fname.size() == 0)
        log.info("Sample file: None, use all samples");
    else
        log.info("Sample file: %s", samp_fname.c_str());

    log.info("Output matrix file: %s", out_fname.c_str());

    int status = FAILED_CALC;

    bcfio::ReadBcf bfid { vcf_fname.c_str() };
    Matrix cov { bfid.n_samples(), bfid.n_samples() };

    if (use_gt) {
        log.info("Relationship matrix: genotype");
        status = compute_genotype_matrix();
    } else if (use_ds) {
        log.info("Relationship matrix: expected alt allele count
        status = compute_eac_matrix();
    } else if (use_both) {
        log.info("Relationship matrix: expected alt allele and haplotype"
                " counts");
        status = compute_geno_and_haplo_matrix();
    } else {
        log.info("Relationship matrix: haplotype");
        status = compute_haplotype_matrix(&log, &bfid, &cov);
    }

    if (status == FAILED_CALC)
        log.error("Computation failed");


    return status;
}
