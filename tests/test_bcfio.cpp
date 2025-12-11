
#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <cstdio>

namespace htslib {
extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}
}

#include <bcfio.h>


char VCF_NAME[] { "build/geno_test_data.vcf" };
char VCFGZ_NAME[] { "build/geno_test_data.vcf.gz" };
char BCF_NAME[] { "build/geno_test_data.bcf" };
size_t K_FOUNDERS = 8;
size_t N_SAMPS = 11;


TEST(TestBcfHeader, ConstructorVcfHdr) {
    htslib::htsFile *fid = htslib::hts_open(VCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());

    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("HD", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, K_FOUNDERS);
    EXPECT_EQ(attr.vl_type, BCF_VL_FIXED);
    EXPECT_EQ(attr.type, BCF_HT_REAL);

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, ConstructorVcfGzHdr) {
    htslib::htsFile *fid = htslib::hts_open(VCFGZ_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());

    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("HD", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, K_FOUNDERS);
    EXPECT_EQ(attr.vl_type, BCF_VL_FIXED);
    EXPECT_EQ(attr.type, BCF_HT_REAL);

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, ConstructorBcfHdr) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("HD", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, K_FOUNDERS);
    EXPECT_EQ(attr.vl_type, BCF_VL_FIXED);
    EXPECT_EQ(attr.type, BCF_HT_REAL);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFmtGt) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("GT", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, 1);
    EXPECT_EQ(attr.vl_type, BCF_VL_FIXED);
    EXPECT_EQ(attr.type, BCF_HT_STR);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFmtGp) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("GP", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, 3);
    EXPECT_EQ(attr.vl_type, BCF_VL_FIXED);
    EXPECT_EQ(attr.type, BCF_HT_REAL);

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, BcfHdrFmtDs) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("DS", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, 1);
    EXPECT_EQ(attr.vl_type, BCF_VL_FIXED);
    EXPECT_EQ(attr.type, BCF_HT_REAL);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFmtErr) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format("DOESNOTEXIST", &attr);
    EXPECT_NE(status, 0);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFilter) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_filter("PASS", &attr);
    EXPECT_EQ(status, 0);

    status = hdr.get_filter("PASSING", &attr);
    EXPECT_NE(status, 0);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrInfoEaf) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_info("EAF", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.type, BCF_HT_REAL);
    EXPECT_EQ(attr.vl_type, BCF_VL_VAR);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrInfoErc) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_info("ERC", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.type, BCF_HT_REAL);
    EXPECT_EQ(attr.vl_type, BCF_VL_VAR);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrInfoErr) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_info("NOTAINFOMEMBER", &attr);
    EXPECT_NE(status, 0);

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, BcfHdrNull) {
    htslib::htsFile *fid = htslib::hts_open("doesnotexist", "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_TRUE(hdr.isnull());
    if (fid) htslib::hts_close(fid);
}



TEST(TestReadBcf, Constructor) {
    bcfio::ReadBcf bcf { VCF_NAME };
    EXPECT_EQ(bcf.n_samples(), N_SAMPS);
    EXPECT_EQ(bcf.k_founders(), K_FOUNDERS);
}

TEST(TestReadBcf, VcfSampNames) {
    bcfio::ReadBcf bcf { VCF_NAME };

    std::unique_ptr<std::string[]> s = bcf.sample_names();

    char samp_name[] = "S01";

    for (int i = 0; i < bcf.n_samples(); i++) {
        snprintf(samp_name, 4, "S%02d", i+1);
        EXPECT_STREQ(s[i].c_str(), samp_name);
    }
}


TEST(TestReadBcf, VcfGzSampNames) {
    bcfio::ReadBcf bcf { VCFGZ_NAME };

    std::unique_ptr<std::string[]> s = bcf.sample_names();

    char samp_name[] = "S01";

    for (int i = 0; i < bcf.n_samples(); i++) {
        snprintf(samp_name, 4, "S%02d", i+1);
        EXPECT_STREQ(s[i].c_str(), samp_name);
    }
}


TEST(TestReadBcf, BcfSampNames) {
    bcfio::ReadBcf bcf { BCF_NAME };

    std::unique_ptr<std::string[]> s = bcf.sample_names();

    char samp_name[] = "S01";

    for (int i = 0; i < bcf.n_samples(); i++) {
        snprintf(samp_name, 4, "S%02d", i+1);
        EXPECT_STREQ(s[i].c_str(), samp_name);
    }
}


// TEST(TestHaplotypeVCFParser, LoadRecord) {
// 
//     HaplotypeVcfParser vcf { VCF_NAME };
// 
//     HaplotypeDataRecord record { vcf.n_samples(), vcf.k_founders() };
// 
//     bool record_loaded { false };
//     record_loaded = vcf.load_record(record);
// 
//     EXPECT_TRUE(record_loaded);
// 
//     EXPECT_EQ(record.chrom(), "chr12");
//     EXPECT_EQ(record.pos(), 788);
//     EXPECT_EQ(record.id(), ".");
//     EXPECT_EQ(record.ref(), 'A');
//     EXPECT_EQ(record.alt(), 'G');
//     EXPECT_EQ(record.qual(), ".");
//     EXPECT_EQ(record.filter(), "PASS");
//     EXPECT_EQ(record.info(), "EAF=0.00228;INFO_SCORE=1;HWE=1;ERC=0.01949;EAC=7.94153;PAF=0.00245;REF_PANEL=0");
//     EXPECT_EQ(record.format(), "GT:GP:DS:HD");
// 
//     EXPECT_EQ(record(0,0), 1.004);
//     EXPECT_EQ(record(0,1), 0);
//     EXPECT_EQ(record(0,2),0.002);
//     EXPECT_EQ(record(0,3),0.001);
//     EXPECT_EQ(record(0,4),0);
//     EXPECT_EQ(record(0,5),0.991);
//     EXPECT_EQ(record(0,6), 0.001);
//     EXPECT_EQ(record(0,7), 0.002);
// 
//     EXPECT_EQ(record(1,0), 0.001);
//     EXPECT_EQ(record(1,1), 0);
//     EXPECT_EQ(record(1,2),0);
//     EXPECT_EQ(record(1,3),0.989);
//     EXPECT_EQ(record(1,4),0.005);
//     EXPECT_EQ(record(1,5),0.005);
//     EXPECT_EQ(record(1,6), 1);
//     EXPECT_EQ(record(1,7), 0);
// 
//     EXPECT_EQ(record(2,0), 0.998);
//     EXPECT_EQ(record(2,1), 0);
//     EXPECT_EQ(record(2,2),0.001);
//     EXPECT_EQ(record(2,3),0);
//     EXPECT_EQ(record(2,4),0);
//     EXPECT_EQ(record(2,5),0);
//     EXPECT_EQ(record(2,6), 1);
//     EXPECT_EQ(record(2,7), 0);
// 
//     EXPECT_EQ(record(3,0), 0.84);
//     EXPECT_EQ(record(3,1), 0);
//     EXPECT_EQ(record(3,2),0);
//     EXPECT_EQ(record(3,3),0);
//     EXPECT_EQ(record(3,4),0);
//     EXPECT_EQ(record(3,5),0);
//     EXPECT_EQ(record(3,6), 1.159);
//     EXPECT_EQ(record(3,7), 0);
// 
//     EXPECT_EQ(record(10,0), 0.592);
//     EXPECT_EQ(record(10,1), 0);
//     EXPECT_EQ(record(10,2),1);
//     EXPECT_EQ(record(10,3),0);
//     EXPECT_EQ(record(10,4),0);
//     EXPECT_EQ(record(10,5),0);
//     EXPECT_EQ(record(10,6),0);
//     EXPECT_EQ(record(10,7),0.407);
// 
//     // load second record
//     record_loaded = vcf.load_record(record);
//     EXPECT_EQ(record.chrom(), "chr12");
//     EXPECT_EQ(record.pos(), 1321);
//     EXPECT_EQ(record.id(), ".");
//     EXPECT_EQ(record.ref(), 'A');
//     EXPECT_EQ(record.alt(), 'C');
//     EXPECT_EQ(record.qual(), ".");
//     EXPECT_EQ(record.filter(), "PASS");
//     EXPECT_EQ(record.info(), "EAF=0.01487;INFO_SCORE=0.17212;HWE=1;ERC=1.33325;EAC=116.998;PAF=0.01127;REF_PANEL=0");
//     EXPECT_EQ(record.format(), "GT:GP:DS:HD");
// 
//     EXPECT_EQ(record(0,0),1.004);
//     EXPECT_EQ(record(0,1), 0);
//     EXPECT_EQ(record(0,2),0.002);
//     EXPECT_EQ(record(0,3),0.001);
//     EXPECT_EQ(record(0,4),0);
//     EXPECT_EQ(record(0,5),0.991);
//     EXPECT_EQ(record(0,6), 0.001);
//     EXPECT_EQ(record(0,7), 0.002);
// 
//     EXPECT_EQ(record(10,0), 0.592);
//     EXPECT_EQ(record(10,1), 0);
//     EXPECT_EQ(record(10,2),1);
//     EXPECT_EQ(record(10,3),0);
//     EXPECT_EQ(record(10,4),0);
//     EXPECT_EQ(record(10,5),0);
//     EXPECT_EQ(record(10,6),0);
//     EXPECT_EQ(record(10,7),0.407);
// 
// }
