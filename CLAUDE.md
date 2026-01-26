# Purpose
- Compute the genetic similarity matrix from genetic data stored in the 
vcf, vcf.gz, and bcf files.
- Use htslib C library to query / read genentic data.
- Write data to a binary file with a header and payload.
    * header data must include all information that is
    required to reproduce the calculation.
    * payload consists of the computed values
- Easy to use command line interface


# External libraries
- htslib, code available on github at: https://github.com/samtools/htslib
- argparse, code available on github at: https://github.com/robert-vogel/argparse


# Code style
- variable, function, class, etc. names:
    * should be descriptive (self-documenting) and short.
    * classes and structs use Upper camel case
    * functions and variables use snake case
    - global constants use all upper case
- keep orthogonal services of the code in distinct modules
- module header files 
    * written to the 'include' directory.
    * API documentation should be on line(s) proceeding the entities they
    describe.
- implementation files in the 'src' directory, make sure to describe the
purpose / design of complex code blocks in the comments
- Unit tests should be written in the 'tests' directory and use Google
test and mock frameworks


# File and directory structure
- 'scratch' directory is not tracked by version control and is for testing
ideas by implementing in small programs
- 'build' directory is the target path for program and unit test builds
and should not be tracked by version control.
- 'include' directory for header files
- 'src' directory for 'main.cpp' file and module implementation files
