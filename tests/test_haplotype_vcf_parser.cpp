
#include "../include/HaplotypeVcfParser.h"
#include "../include/utils.h"
#include <gtest/gtest.h>



char VCF_NAME[] { "../tests/test.vcf" };


TEST(TestHaplotypeVCFParser, Constructor) {

    HaplotypeVcfParser vcf { VCF_NAME };
    
    EXPECT_EQ(vcf.n_samples(), 11);
    EXPECT_EQ(vcf.k_founders(), 8);
}


TEST(TestHaplotypeVCFParser, LoadRecord) {

    HaplotypeVcfParser vcf { VCF_NAME };

    HaplotypeDataRecord record { vcf.n_samples(), vcf.k_founders() };

    bool record_loaded { false };
    record_loaded = vcf.load_record(record);

    EXPECT_TRUE(record_loaded);

    EXPECT_EQ(record.chrom(), "chr12");
    EXPECT_EQ(record.pos(), 788);
    EXPECT_EQ(record.id(), ".");
    EXPECT_EQ(record.ref(), 'A');
    EXPECT_EQ(record.alt(), 'G');
    EXPECT_EQ(record.qual(), ".");
    EXPECT_EQ(record.filter(), "PASS");
    EXPECT_EQ(record.info(), "EAF=0.00228;INFO_SCORE=1;HWE=1;ERC=0.01949;EAC=7.94153;PAF=0.00245;REF_PANEL=0");
    EXPECT_EQ(record.format(), "GT:GP:DS:HD");

    EXPECT_EQ(record(0,0), 1.004);
    EXPECT_EQ(record(0,1), 0);
    EXPECT_EQ(record(0,2),0.002);
    EXPECT_EQ(record(0,3),0.001);
    EXPECT_EQ(record(0,4),0);
    EXPECT_EQ(record(0,5),0.991);
    EXPECT_EQ(record(0,6), 0.001);
    EXPECT_EQ(record(0,7), 0.002);

    EXPECT_EQ(record(1,0), 0.001);
    EXPECT_EQ(record(1,1), 0);
    EXPECT_EQ(record(1,2),0);
    EXPECT_EQ(record(1,3),0.989);
    EXPECT_EQ(record(1,4),0.005);
    EXPECT_EQ(record(1,5),0.005);
    EXPECT_EQ(record(1,6), 1);
    EXPECT_EQ(record(1,7), 0);

    EXPECT_EQ(record(2,0), 0.998);
    EXPECT_EQ(record(2,1), 0);
    EXPECT_EQ(record(2,2),0.001);
    EXPECT_EQ(record(2,3),0);
    EXPECT_EQ(record(2,4),0);
    EXPECT_EQ(record(2,5),0);
    EXPECT_EQ(record(2,6), 1);
    EXPECT_EQ(record(2,7), 0);

    EXPECT_EQ(record(3,0), 0.84);
    EXPECT_EQ(record(3,1), 0);
    EXPECT_EQ(record(3,2),0);
    EXPECT_EQ(record(3,3),0);
    EXPECT_EQ(record(3,4),0);
    EXPECT_EQ(record(3,5),0);
    EXPECT_EQ(record(3,6), 1.159);
    EXPECT_EQ(record(3,7), 0);

    EXPECT_EQ(record(10,0), 0.592);
    EXPECT_EQ(record(10,1), 0);
    EXPECT_EQ(record(10,2),1);
    EXPECT_EQ(record(10,3),0);
    EXPECT_EQ(record(10,4),0);
    EXPECT_EQ(record(10,5),0);
    EXPECT_EQ(record(10,6),0);
    EXPECT_EQ(record(10,7),0.407);

    // load second record
    record_loaded = vcf.load_record(record);
    EXPECT_EQ(record.chrom(), "chr12");
    EXPECT_EQ(record.pos(), 1321);
    EXPECT_EQ(record.id(), ".");
    EXPECT_EQ(record.ref(), 'A');
    EXPECT_EQ(record.alt(), 'C');
    EXPECT_EQ(record.qual(), ".");
    EXPECT_EQ(record.filter(), "PASS");
    EXPECT_EQ(record.info(), "EAF=0.01487;INFO_SCORE=0.17212;HWE=1;ERC=1.33325;EAC=116.998;PAF=0.01127;REF_PANEL=0");
    EXPECT_EQ(record.format(), "GT:GP:DS:HD");

    EXPECT_EQ(record(0,0),1.004);
    EXPECT_EQ(record(0,1), 0);
    EXPECT_EQ(record(0,2),0.002);
    EXPECT_EQ(record(0,3),0.001);
    EXPECT_EQ(record(0,4),0);
    EXPECT_EQ(record(0,5),0.991);
    EXPECT_EQ(record(0,6), 0.001);
    EXPECT_EQ(record(0,7), 0.002);

    EXPECT_EQ(record(10,0), 0.592);
    EXPECT_EQ(record(10,1), 0);
    EXPECT_EQ(record(10,2),1);
    EXPECT_EQ(record(10,3),0);
    EXPECT_EQ(record(10,4),0);
    EXPECT_EQ(record(10,5),0);
    EXPECT_EQ(record(10,6),0);
    EXPECT_EQ(record(10,7),0.407);

}
