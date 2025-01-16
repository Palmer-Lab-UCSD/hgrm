
# By: Robert Vogel
# Affiliation: Palmer Lab at UCSD
# 2025-01-11
# 
# Acknowledgment
# This file was originally written by Robert Vogel.  Claude,
# the AI assistant by Anthropic reviewed this file and provided
# recommendations.

cmake_minimum_required(VERSION 3.31)
project(hgrm VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions")

# Begin Claude recommendation
# include_directories(${PROJECT_SOURCE_DIR}/include)
# 
# add_library(hgrm_lib
#     src/HaplotypeDataRecord.cpp
#     src/Matrix.cpp
#     src/HaplotypeVcfParser.cpp
# )
# 
# add_executable(hgrm src/main.cpp)
# target_link_libraries(hgrm PRIVATE hgrm_lib)

# end Claude recommendation


find_package(GTest REQUIRED)

add_library(matrix_lib src/Matrix.cpp)
target_include_directories(matrix_lib PUBLIC include)

add_library(parse_lib src/HaplotypeDataRecord.cpp src/HaplotypeVcfParser.cpp)
target_include_directories(parse_lib PUBLIC include)



# Testing configuration
enable_testing()

add_executable(
    test_matrix
    tests/test_matrix.cpp
)

target_link_libraries(
    test_matrix
    PRIVATE
    matrix_lib
    GTest::gtest_main
)

add_executable(
    test_haplotype_data_record
    tests/test_haplotype_data_record.cpp
)

add_executable(
    test_haplotype_vcf_parser
    tests/test_haplotype_vcf_parser.cpp
)

target_link_libraries(
    test_haplotype_data_record
    PRIVATE
    matrix_lib
    parse_lib
    GTest::gtest_main
)


target_link_libraries(
    test_haplotype_vcf_parser
    PRIVATE
    matrix_lib
    parse_lib
    GTest::gtest_main
)


# add_executable(
#     hgrm
#     src/main.cpp
# )
# 
# target_link_libraries(
#     hgrm
#     PRIVATE
#     matrix_lib
#     parse_lib
# )


include(GoogleTest)
gtest_discover_tests(test_matrix)
gtest_discover_tests(test_haplotype_data_record)
gtest_discover_tests(test_haplotype_vcf_parser)

