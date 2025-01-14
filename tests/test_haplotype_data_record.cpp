

#include "../include/HaplotypeVcfParser.h"
#include <gtest/gtest.h>




TEST(TestConstructorAssignment, Constructor) {

    size_t num_vcf_columns { 13 };
    size_t num_samps { 4 };
    size_t num_founders { 3 };
    static constexpr int num_records { 3 };


    std::string vcf_record { "chr12 1 . A T Q1 F1 INFO1 GT:AB:HD 0/0:2:1,0,1 0/0:1:0,2,0 1/0:0:0,0,2 1/1:1:2,0,0\n" };


    // Test first record
    HaplotypeDataRecord hap_record { vcf_record, num_vcf_columns };
    
    std::array<size_t, 2> dims { hap_record.dims() };
    EXPECT_EQ(dims[0], num_founders);
    EXPECT_EQ(dims[1], num_samps);
    
    EXPECT_EQ(hap_record(0,0), 1);
    EXPECT_EQ(hap_record(0,1), 0);
    EXPECT_EQ(hap_record(0,2), 0);
    EXPECT_EQ(hap_record(0,3), 2);
    
    EXPECT_EQ(hap_record(1,0), 0);
    EXPECT_EQ(hap_record(1,1), 2);
    EXPECT_EQ(hap_record(1,2), 0);
    EXPECT_EQ(hap_record(1,3), 0);

    EXPECT_EQ(hap_record(2,0), 1);
    EXPECT_EQ(hap_record(2,1), 0);
    EXPECT_EQ(hap_record(2,2), 2);
    EXPECT_EQ(hap_record(2,3), 0);

    // Test second record
     
    vcf_record = "chr12 2 . G T Q2 F2 INFO2 GT:AB:HD 0/0:2:1,0,1 0/1:1:0,2,0 1/0:0:0,0,2 0/1:1:0,0,2\n";
    HaplotypeDataRecord hap_record2 { vcf_record, num_vcf_columns };

    EXPECT_EQ(hap_record2(0,0), 1);
    EXPECT_EQ(hap_record2(0,1), 0);
    EXPECT_EQ(hap_record2(0,2), 0);
    EXPECT_EQ(hap_record2(0,3), 0);
    
    EXPECT_EQ(hap_record2(1,0), 0);
    EXPECT_EQ(hap_record2(1,1), 2);
    EXPECT_EQ(hap_record2(1,2), 0);
    EXPECT_EQ(hap_record2(1,3), 0);

    EXPECT_EQ(hap_record2(2,0), 1);
    EXPECT_EQ(hap_record2(2,1), 0);
    EXPECT_EQ(hap_record2(2,2), 2);
    EXPECT_EQ(hap_record2(2,3), 2);
    
}


TEST(TestConstructorAssignment, DefaultConstructor) {
    HaplotypeDataRecord data {};

    EXPECT_EQ(data.chrom(), "");
    EXPECT_EQ(data.pos(), -1);
    EXPECT_EQ(data.id(), "");
    EXPECT_EQ(data.ref(), '\0');
    EXPECT_EQ(data.alt(), '\0');
    EXPECT_EQ(data.qual(), "");
    EXPECT_EQ(data.filter(), "");
    EXPECT_EQ(data.info(), "");
    EXPECT_EQ(data.format(), "");
}


TEST(TestConstructorAssignment, CopyConstructor) {
    
    size_t num_vcf_columns { 13 };


    std::string vcf_record { "chr12 1 . A T Q1 F1 INFO1 GT:AB:HD 0/0:2:1,0,1 0/0:1:0,2,0 1/0:0:0,0,2 1/1:1:2,0,0\n" };

    HaplotypeDataRecord origin { vcf_record, num_vcf_columns };
    HaplotypeDataRecord replicate { origin };

    // test that the objects are indeed different
    EXPECT_NE(&origin, &replicate);
    
    std::array<size_t, 2> o_dims { origin.dims() };
    std::array<size_t, 2> rep_dims { replicate.dims() };

    EXPECT_EQ(o_dims[0], rep_dims[0]);
    EXPECT_EQ(o_dims[1], rep_dims[1]);
    

    // make sure that copy constructor correct
    EXPECT_EQ(origin(0,0), replicate(0,0));
    EXPECT_EQ(origin(0,1), replicate(0,1));
    EXPECT_EQ(origin(0,2), replicate(0,2));
    EXPECT_EQ(origin(0,3), replicate(0,3));
    

    EXPECT_EQ(origin(1,0), replicate(1,0));
    EXPECT_EQ(origin(1,1), replicate(1,1));
    EXPECT_EQ(origin(1,2), replicate(1,2));
    EXPECT_EQ(origin(1,3), replicate(1,3));


    EXPECT_EQ(origin(2,0), replicate(2,0));
    EXPECT_EQ(origin(2,1), replicate(2,1));
    EXPECT_EQ(origin(2,2), replicate(2,2));
    EXPECT_EQ(origin(2,3), replicate(2,3));
    
}

//    std::array<std::string,3> meta { "##version=4.2\n",
//                                    "##ID=<test>\n",
//                                    "##FORMAT=<GT>\n" };
//
//    std::string h { "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4\n" };
