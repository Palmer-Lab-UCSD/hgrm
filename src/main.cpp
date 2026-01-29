// Compute the genomic relationship matrix using haplotypes
//
// Palmer Lab at UCSD
//
// This program performs a single-pass computation of the genetic relationship
// matrix (GRM, GR matrix).  GR matrices may be constructed using alt allele
// counts, expected alt allele counts, expected haplotype counts, or a
// combination of both expected alt allele and haplotype counts.  
//
#include <argparse.h>
#include <optional>
#include <string>

#include <logger.h>
#include <calc.h>
#include <grm.h>


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
        "grm: Genetic Relationship Matrix",
        "This program provides tools for computing the genetic relationship"
        " matrix (GRM) and the leave-one-chromosome-out (LOCO) matrices for"
        " linear mixed effect based association studies.  The GRM may be"
        " computed using called genotypes, expected alternative allele counts,"
        " expected haplotype counts, or both expected alternative allele"
        " and haplotype counts.  By default the expected alternative allele"
        " counts are used."
    };

    argparse::CmdDef *contig_cmd = parser.add_cmd("contig");

    contig_cmd->add_arg("-o", 
            argparse::ArgType::STRING,
            "the path and filename that the resulting haplotype genetic"
            " relationship matrix is printed.");

    contig_cmd->add_arg("--sample_names",
            argparse::ArgType::STRING,
            "The path and name of the file containing sample names to be"
            " included in computing the relationship matrix.  The file must"
            " include a single sample filename, and if necessary file system"
            " path, per line.");

    contig_cmd->add_arg("--gt",
            argparse::ArgType::BOOLEAN,
            "Use sample genotypes to compute the relationship matrix");

    contig_cmd->add_arg("--ehc",
            argparse::ArgType::BOOLEAN,
            "Use sample expected haplotype count to compute the the genetic"
            " relationship matrix.");

    contig_cmd->add_arg("-b",
            argparse::ArgType::BOOLEAN,
            "Use both the expected alternative allele and haplotype counts to"
            " compute the genetic relationship matrix");

    contig_cmd->add_arg("bcf",
            argparse::ArgType::STRING, 
            "The path and filename of the genetic data to compute the GRM. The"
            " data may be in any of the htslib supported formats, i.e. vcf,"
            " vcf.gz, or bcf.");


    argparse::CmdDef *loco_cmd = parser.add_cmd("loco");
    loco_cmd->add_arg("filename",
            argparse::ArgType::STRING,
            "Name, and path, of file that stores the name and paths of matrix"
            " files used to compute leave-one-chromosome-out (LOCO) relationship"
            " matrix.");



    Logger log {};
    int status = FAILED_CALC;
    argparse::ArgStatus arg_status = parser.parse_args(argc, argv);

    // PARSE ARGUMENTS
    if (arg_status == argparse::ArgStatus::HELP)
        return 0;

    if (arg_status != argparse::ArgStatus::SUCCESS) {
        log.error("Error: couldn't parse command line args, exiting\n");
        exit(EXIT_FAILURE);
    }

    // EXTRACT ARGS
    if (parser.is_sub_cmd("contig")) {

        std::optional<std::string> tmp_str {};
        if((tmp_str = parser.get<std::string>("bcf")) == std::nullopt) {
            log.error("Error retrieving vcf name");
            exit(EXIT_FAILURE);
        }
        std::string bcf_fname { tmp_str.value() };

        if ((tmp_str = parser.get<std::string>("o")) == std::nullopt) {
            log.error("Error retrieving output name");
            exit(EXIT_FAILURE);
        }
        std::string out_fname { tmp_str.value() };

        if (out_fname.size() == 0)
            out_fname = bcf_fname + ".mat";

        if ((tmp_str = parser.get<std::string>("sample_names")) == std::nullopt) {
            log.error("Error retrieving sample_names file.\n");
            exit(EXIT_FAILURE);
        }

        std::string samp_fname { tmp_str.value() };

        
        std::optional<bool> tmp_bool {};
        if ((tmp_bool = parser.get<bool>("gt")) == std::nullopt) {
            log.error("Error retrieving relationship matrix type.\n");
            exit(EXIT_FAILURE);
        }
        bool use_gt { tmp_bool.value() };

        if ((tmp_bool = parser.get<bool>("ehc")) == std::nullopt) {
            log.error("Error retrieving relationship matrix type.\n");
            exit(EXIT_FAILURE);
        }
        bool use_ehc { tmp_bool.value() };

        if ((tmp_bool = parser.get<bool>("b")) == std::nullopt) {
            log.error("Error retrieving relationship matrix type.\n");
            exit(EXIT_FAILURE);
        }
        bool use_both { tmp_bool.value() };



        if ((use_gt && use_both) || (use_gt && use_ehc) || (use_both && use_ehc)) {
            log.error("user must specify either use_gt, use_both, use_ds,"
                    " or omit both options to compute the haplotype based"
                    " relationship matrix.");
            exit(EXIT_FAILURE);
        }
    
        log.info("BCF/VCF file name: %s", bcf_fname.c_str());


        bcfio::ReadBcf bfid { bcf_fname.c_str() };

        int bstatus = 0;
        if (samp_fname.size() == 0)
            log.info("Sample file: None, use all samples");
        else if ((bstatus = bfid.set_samples(samp_fname.c_str())) == 0)
            log.info("Sample file: %s", samp_fname.c_str());
        else if (bstatus < 0) {
            log.error("Subsetting by sample file, %s, resulted in error", 
                    samp_fname.c_str());
            return -1;
        } else if (bstatus > 0) {
            log.error("One or more samples specified in sample file, %s,"
                    " do not %s", 
                    samp_fname.c_str(), 
                    bcf_fname.c_str());
            return -1;
        }

        log.info("Output matrix file: %s", out_fname.c_str());

        grm.Grm cov { bfid.n_samples(), bfid.n_samples() };

        if (use_gt) {
            log.info("Relationship matrix: genotype");
            status = compute_genotype_matrix();
        } else if (use_ehc) {
            log.info("Relationship matrix: expected haplotype count");
            status = compute_ehc_matrix(&log, &bfid, &cov);
        } else if (use_both) {
            log.info("Relationship matrix: expected alt allele and haplotype"
                    " counts");
            status = compute_eac_and_ehc_matrix();
        } else {
            log.info("Relationship matrix: expected alternative allele counts");
            status = compute_eac_matrix();
        }

        if (status == FAILED_CALC)
            log.error("Computation failed");

        log.info("Writing to file");

        cov.write(out_fname);
    }

    if (parser.is_sub_cmd("loco"))
        printf("loco selected\n");


    return status;
}
