# `grm` a tool for computing genetic relationship matrices

 🏗️  **Under construction** 🏗️


## Table of Contents

1. [About](#about)
1. [Subprograms: Compute GRM and leave-one-chromosome-out](#subprog)
2. [Command line user interface](#cli)
3. [Installation and requirements](#install)
4. [References](#refs)

# About <a name="about"></a>

The genetic relationship matrix (GRM) describes the genetic relationship 
between pairs of samples.  Its computation is dependent on the random effects
defined by the linear mixed model mapping genetic features to phenotype.  For
example, suppose that we are interested in accounting for polygenic SNP effects
using measured genotypes.  Let the number of samples be $N$, the number of
loci that contribute to polygenic effects $M$, and the quantiative phenotypes
of $N$ samples $Y\in\mathbb{R}^{N\times 1}$. Under the LMM [[1]](#refs)

$$
Y = x_j\beta_j + \mathbf{Z}_j U_j + \epsilon.\\
$$

with the fixed effect at locus $j$ being the alternative allele count, denoted
$x_j\in \{0,1,2\}^{N\times 1}$, and fixed effect size $\beta_j$.  The
random polygenic and environmental effects have properties

$$
\begin{align}
U_j &\sim \mathcal{N}\left(0, \sigma_g^2\; \mathbf{I}_{M\times M}\right)\\
\epsilon &\sim \mathcal{N}\left(0, \sigma_e^2\; \mathbf{I}_{N\times N}\right)
\end{align}
$$

with $\mathbf{Z}_j\in\{0,1,2\}^{N\times M}$ being the matrix of alt allele
counts of the $N$ samples and the set of $M$ markers in which locus $j$ is
not a member.

Under this model the sample phenotype covariance matrix
decomposes into genetic and environmental terms

$$
\text{cov}(Y) = \overbrace{\sigma^2_g\mathbf{Z}\mathbf{Z}^T}^{\text{genetic}}   
+
\underbrace{\sigma^2_e\mathbf{I}_{N\times N}}_{\text{Environment}}
$$

where the genetic component of the phenotype covariance tells us how
to compute the GRM, i.e. $\mathbf{Z}\mathbf{Z}^T$.  

In general, the alt allele count polygenic random effects are not the only 
genetic effects that we may account for.  This program includes GRMs 
computed from genotypes, the expected alt allele count under the imputation
model, the expected haplotype count, and a combination of the expected alt 
allele count with the expected haplotype count.


## Genetic relationship matrices

The genetic relationship matrices modeling distinct genetic random effects
are derived as outlined above.  Here we enumerate the GRM for each type
of random effect consider.


### SNP GRM

The SNP GRM is presented in the [about](#about) section.  Let $A_\text{SNP}$ be the
GRM computed by polygenic SNP effects, then 

$$
G = ZZ^T
$$

with $Z$ being the $N\times M$ matrix of alt allele counts.


### Expected alt allele count GRM

In many cases, as is the case in the Palmer Lab, SNP are imputed.  If the
imputation method estimates the genotype probabilities of each locus of each
sample as out method of choice, STITCH [[2]](#refs), then we are able to
compute the expected alt allele count

$$
\begin{align}
\mathbb{E}[X_{ji} | p_{ji}] = \sum_{x_{ji} = 0}^3 x_{ji}\;\mathbb{P}(X_{ji})
\end{align}
$$

where $p

### Expected haplotype count GRM

## Subprograms: Compute GRM and leave-one-chromosome-out <a name="#subprog"></a>

The `grm` program consists of two subprograms: 

* `grm contig`: the computation of the GRM of a named contig
* `grm loco`: the aggregation of contig GRMs into a 
    "leave-one-chromosome-out" matrix (denoted the "loco" matrix).  

each of which have there own options enumerated in the next section.4=



## Installation and requirements <a name="install"></a>

The program is only available as source from this repository and requires

* `GNU make` 
* `htslib` https://github.com/samtools/htslib
* `argparse` https://github.com/robert-vogel/argparse
* `clang` or `gcc` C++17 compiler


## Genetic Relationship Matrix types



## Compute the genetic relationship matrix

The GRM calculation requires the SNPs or haplotypes to jbe in the bcf family
of file formats, i.e. vcf, vcf.gz, or bcf.  By default, the GRM is computed
using the expected haplotype counts with FORMAT ID = "HD".

```
grm chrm <chrm_id> <path/to/my/snps.bcf>
```

will produce a binary `.mat` file that stores the GRM and relavent meta data.
Other options include




## Compute LOCO matrices

```
grm loco path/to/file/with/grm_filename_and_path_per_line
```





## Contributing

I am using [GoogleTest](https://google.github.io/googletest/) framework
for organizing tests.  If you contribute, please make tests for your
contributions.  To run tests, `build` directory and build the project
```
make check


## Acknowledgement

Code design and original version completed by Robert Vogel,
reviewed by Claude Sonnet, the AI assistant from Anthropic
(Jan 2025), with minor recommendations incorporated.

## References <a name="refs"></a>

[1] [Yang et al. Nature Genetics 42, 565-569 (2010)](https://www.nature.com/articles/ng.608)

[1] [Kang et al. Genetics 178: 1709-1723 (2008)](https://academic.oup.com/genetics/article/178/3/1709/6061473)

[2] [Kang et al. Nature Genetics 42 348-354 (2010)](https://www.nature.com/articles/ng.548)


