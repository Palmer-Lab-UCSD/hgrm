# Compute the genetic relationship matrix using expected haplotype counts



The genetic relationship matrix (GRM) is the covariance between samples over
the measured genetic markers.  It's utility is in accounting for 
relatedness among samples, as random effects in a Linear Mixed Effects
Model, when performing a Genome Wide Association Study (GWAS) [1] and computing
the heritability of complex traits [2].  The GRM is traditionally computed using
SNP genotypes, but here we are interested in haplotype based covariance.

This software is aimed at computing N sample by N sample covariance 
matrix over the expected haplotype counts and markers.  Let
$x_{mpi}\in[0,2]$ be the expected count of marker $m$ for 
founder $p$ and sample $i$.  The algorithm we will use to compute the
covariance is an online, or single-pass,
algorithm [3] modified for our expected haplotype count encoding.  Let 
$m$ indicated we are computing the update of the $m^{th}$ marker given
computations of the mean \$\bar{x}\_{(m-1)pi}\$ and covariance between 
samples \$i\$ and \$j\$, and \$C\_{ij}^{(m-1)}\$ of the preceeding $m-1$ markers.
The recursive mean and covariance update
equations are as follows

$$
\bar{x}\_{mpi} = \bar{x}\_{(m-1)pi} + \frac{\delta_{mpi}}{m}
$$

where

$$
\delta_{mpi} = x_{mpi} - \bar{x}_{(m-1)pi}.
$$

The update to the covariance is

$$
C_{ij}^{(m)} = \frac{M-2}{M-1} \\; C_{ij}^{(M-1)} 
    + \frac{1}{M}\sum_{p=1}^K \delta_{Mpi}\delta_{Mpj}.
$$


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

[1] [Kang et al. Nature Genetics 42 348-354 (2010)](https://www.nature.com/articles/ng.548)

[2] [Yang et al. Nature Genetics 42, 565-569 (2010)](https://www.nature.com/articles/ng.608)

[3] [Wikipedia, online covariance algorithm](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online)
