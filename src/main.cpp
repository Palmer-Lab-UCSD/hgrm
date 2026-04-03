// Haplotype wide association analyses
//
// Palmer Lab at UCSD
//
// This program provides tools to associate ancestor haplotype at
// any locus to the phenotype of interests.  We assume that the
// expected ancestor haplotype count is provided in the vcf, vcf.gz,
// or bcf file formats.  The analyses conducted are:
//
// * computation of the haplotype genetic relationship matrix
// * estimation of the genetic and environment variances 
// * estimation of haplotype heritability
// * computation of log odds (lod) score at each locus
// * computation of the ancestor Best Linear Unbiased Estimator
//  (BLUPs) with phenotype at each locus.
//
#include <argparse.h>
#include <optional>
#include <string>

#include <logger.hpp>
#include <grm.hpp>


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
        "hwas: Haplotype-wide association analyses",
        " This program provides tools for computing statistical"
        " associations between ancestor founder haplotype at any"
        " locus i with a quantitative phenotype under a linear"
        " mixed effects model (LMM)."
    };

    argparse::CmdDef *grm_cmd = parser.add_cmd("grm");

    grm_cmd->add_arg("-o", 
            argparse::ArgType::STRING,
            "the path and filename that the resulting haplotype"
            " genetic relationship matrix is stored.");

    grm_cmd->add_arg("--samples",
            argparse::ArgType::STRING,
            "The path and name of the file containing sample names"
            " to be included in computing the relationship matrix."
            " The file must include a single sample filename, and if"
            " necessary file system path, per line.");

    grm_cmd->add_arg("--eac",
            argparse::ArgType::BOOLEAN,
            "Use sample expected alternative allele counts to compute"
            " the relationship matrix.");

    grm_cmd->add_arg("--ehc",
            argparse::ArgType::BOOLEAN,
            "Use sample expected haplotype count to compute the the"
            " genetic relationship matrix.");

    grm_cmd->add_arg("-b",
            argparse::ArgType::BOOLEAN,
            "Use both the expected alternative allele and haplotype"
            " counts to compute the genetic relationship matrix");

    grm_cmd->add_arg("bcf",
            argparse::ArgType::STRING, 
            "The path and filename of the genetic data to compute the"
            " GRM. The data may be in any of the htslib supported"
            " formats, i.e. vcf, vcf.gz, or bcf.");


    argparse::CmdDef *loco_cmd = parser.add_cmd("loco");
    loco_cmd->add_arg("filename",
            argparse::ArgType::STRING,
            "Name, and path, of file that stores the name and paths of"
            " matrix files used to compute leave-one-chromosome-out"
            " (LOCO) relationship matrix.");

    argparse::CmdDef *lmm_cmd = parser.add_cmd("lmm");
    lmm_cmd->add_arg("bcf",
            argparse::ArgType::STRING,
            "Name, and path, of vcf, vcf.gz, or bcf file with the"
            " genotype data.");
    lmm_cmd->add_arg("--samples",
            argparse::ArgType::STRING,
            "Name, and path, to the file with sample id's.  This"
            " is used to select data from that stored in the vcf,"
            " vcf.gz, or bcf file.");


    argparse::CmdDef *eqtl_cmd = parser.add_cmd("eqtl");
    eqtl_cmd->add_arg("bcf",
            argparse::ArgType::STRING,
            "Name, and path, of vcf, vcf.gz, or bcf file with the"
            " genotype data.");


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

    // PARSE ARGS FOR RESPECTIVE SUBPROGRAMS AND RUN
    //
    // Compute the GRM for the specified contig
    if (parser.is_sub_cmd("grm")) {

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
            log.error("user must specify either use_gt, use_both, use_ehc,"
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
            exit(EXIT_FAILURE);
        } else if (bstatus > 0) {
            log.error("One or more samples specified in sample file, %s,"
                    " do not %s", 
                    samp_fname.c_str(), 
                    bcf_fname.c_str());
            exit(EXIT_FAILURE);
        }

        log.info("Output matrix file: %s", out_fname.c_str());

        grm::Grm grmatrix { bfid.n_samples() };

        if (use_gt) {
            log.info("Relationship matrix: genotype");
            status = compute_genotype_matrix();
        } else if (use_ehc) {
            log.info("Relationship matrix: expected haplotype count");
            status = compute_ehc_matrix(&log, &bfid, &grmatrix);
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

        grmatrix.write(out_fname);
    } else if (parser.is_sub_cmd("loco")) {
        // Compute the leave-one-chromosome-out matrix given a set of 
        // matricies.
        //
        printf("loco selected\n");
    } else if (parser.is_sub_cmd("lmm")) {
        printf("association statistics selected\n");
    } else if (parser.is_sub_cmd("eqtl")
        printf("eqtl statistics selected\n");


    return status;
}
