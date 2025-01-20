# 🏗️ Being build 🏗️

# Compute the genetic relationship matrix using expected haplotype counts



The genetic relationship matrix (GRM) is the covariance between samples over
the measured genetic markers.  It's utility is in accounting for 
relatedness among samples, as random effects in a Linear Mixed Effects
Model, when performing a Genome Wide Association Study (GWAS) [1,2] and computing
the heritability of complex traits [3].  The GRM is traditionally computed using
SNP genotypes, but here we are interested in haplotype based covariance.


## Running the software

The program is ran by supplying the path and filename of a VCF.  Note, as 
of now this software does not support data streams with UNIX pipe or
bgzip, gzip, etc. compression.  The output is printed to standard out.

```
hgrm path/to/my_vcf > grm
```


## Installation and availability

The program is only available as source from this repository and requires

* `cmake` (>= 3.31.4)
* `make` 
* `clang` or `gcc` C++17 compiler

To install navigate to the top level directory of this repository,
and make the `build` directory
```bash
mkdir build
```
then use `cmake` to generate `make` files, etc.,
```bash
cmake -S . -B build/
```
followed by navigating to the build directory and running GNU `make`
```bash
make
```
In the build directory you should now have a binary file called `hgrm`,
this is the executable program.  Move it to a directory in your
shell's search path.


## Contributing

I am using [GoogleTest](https://google.github.io/googletest/) framework
for organizing tests.  If you contribute, please make tests for your
contributions.  To run tests, `build` directory and build the project
```
cmake -S ../ -B .
make
```
Then use `cmake`'s utility
```
ctest
```

## Acknowledgement

Code design and original version completed by Robert Vogel,
reviewed by Claude Sonnet, the AI assistant from Anthropic
(Jan 2025), with minor recommendations incorporated.

## References

[1] [Kang et al. Genetics 178: 1709-1723 (2008)](https://academic.oup.com/genetics/article/178/3/1709/6061473)

[2] [Kang et al. Nature Genetics 42 348-354 (2010)](https://www.nature.com/articles/ng.548)

[3] [Yang et al. Nature Genetics 42, 565-569 (2010)](https://www.nature.com/articles/ng.608)

