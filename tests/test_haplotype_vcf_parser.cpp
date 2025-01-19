
#include "../include/HaplotypeVcfParser.h"
#include "../include/utils.h"
#include <gtest/gtest.h>



char VCF_NAME[] { "../tests/test.vcf" };


TEST(TestHaplotypeVCFParser, Constructor) {

    HaplotypeVcfParser vcf { VCF_NAME };
    
    EXPECT_EQ(vcf.n_samples(), 11);
    EXPECT_EQ(vcf.k_founders(), 8);

}
