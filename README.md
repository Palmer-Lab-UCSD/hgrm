# `grm` a tool for computing genetic relationship matrices

 🏗️  **Under construction** 🏗️


## Table of Contents

1. [About](#about)
2. [Genetic relationship matrices](#grm)
3. [Command line interface](#cli)
4. [Installation and requirements](#install)
5. [`.grm` file format](#grmspec)
4. [Contributing](#contributing)
4. [A.I. Acknowledgement](#ai)
5. [References](#refs)

## About <a name="about"></a>

The genetic relationship matrix (GRM) describes the genetic relationship 
between pairs of samples.  Its computation is dependent on the random effects
defined by the linear mixed model mapping genetic features to phenotype.  For
example, suppose that we are interested in accounting for polygenic SNP effects
using measured genotypes.  Let the number of samples be $N$, the set of loci of
contributing to polygenic effects be $\Omega$ with $|\Omega|= M$, and the
quantitative phenotypes of $N$ samples $Y\in\mathbb{R}^{N\times 1}$. Note that
locus $j$ is not a member of the set of loci $\Omega$. Under the
LMM [[1]](#refs)

$$
Y = x_j\beta_j + \mathbf{Z} U + \epsilon.\\
$$

with the fixed effect at locus $j$ being the alternative allele count, denoted
$x_j\in \{0,1,2\}^{N\times 1}$, and fixed effect size $\beta_j$.  The
random polygenic and environmental effects have properties

$$
\begin{align}
U &\sim \mathcal{N}\left(0, \sigma_g^2 \\: \mathbf{I}_{M\times M}\right)\\
\epsilon &\sim \mathcal{N}\left(0, \sigma_e^2 \\: \mathbf{I}_{N\times N}\right)
\end{align}
$$

with $\mathbf{Z}\in\{0,1,2\}^{N\times M}$ being the matrix of alt allele
counts of the $N$ samples and and $M$ loci in the set $\Omega$,
not a member.

Under this model the sample phenotype covariance matrix
decomposes into genetic and environmental terms

$$
\text{cov}(Y) = \overbrace{
    \sigma^2_g \\; \mathbf{Z}\mathbf{Z}^T
}^{\text{genetic}} 
+
\underbrace{
    \sigma^2_e \\; \mathbf{I}_{N\times N}
}_{\text{environment}}
$$

where the genetic component of the phenotype covariance tells us how
to compute the GRM, i.e. $\mathbf{Z}\mathbf{Z}^T$.  

The alt allele count polygenic random effects are one example of genetic
effects.  We do not need to limit ourselves to this model, and instead
account for any measurable genetic signals.  This program uses several 
signals to compute a GRM:

* genotypes i.e. the alternative allele count,
* the expected alt allele count under a probabilistic model, useful for
imputed genetic signals,
* ancestral haplotypes, i.e. a K dimensional vector of expected haplotype 
counts,
* or some combination of the aforementioned signals.

This program provides an means to compute the GRM of genetic signals in
general.


## Genetic relationship matrices <a name="grm"></a>

The genetic relationship matrices that this program computes are as follows:
imputed SNP genotypes considered above, the expected alternative allele 
counts under a probabilistic model, the expected ancestral haplotype counts,
and the combined expected alternative allele and expected haplotype counts.
The subsections the follow define the model and the calculation for the
aforementioned GRMs.

### SNP GRM

The SNP GRM is presented in the [about](#about) section.  Let $A_\text{SNP}$ be the
GRM computed by polygenic SNP effects, then 

$$
A_\text{SNP} = \mathbf{Z}\mathbf{Z}^T
$$

with $\mathbf{Z}$ being the $N\times M$ matrix of alt allele counts.


### Expected alternative allele count GRM

In many cases, as is the case in the Palmer Lab, SNPs are ***imputed***.  If
the imputation method estimates the genotype probabilities at each locus of
each sample, as does our method of choice STITCH [[2]](#refs), then it may 
be more appropriate consider the expected alternative allele counts (EAC)
under the imputation model rather than the imputed genotype calls.

Let $Z_{im}\in\{0, 1, 2\}$ be an element of the design matrix of random
polygenic effects for sample $i$ at locus $m\in\Omega$.  We consider $Z_{im}$
to be random as the genotypes are imputed under a probabilistic model. This
makes the expected value $\mathbb{E}[Z_{im} |\mathcal{O}]$ of alternative
allele counts given the experimentally observed reads the an appropriate
genetic signal to consider.  Let's call the matrix of expected alternative
allele counts $\mathbf{C} := \mathbb{E}[\mathbf{Z} |\mathcal{O}]$.  Given
this, the GRM over expected alternative
allele counts $\mathbf{A}_\text{EAC}$ is

$$
\mathbf{A}_\text{EAC} =\mathbf{C}\mathbf{C}^T
$$

### Expected haplotype count GRM

The expected haplotype count GRM $\mathbf{A}_\text{EHC}$ needs
motivation.  The reason is that at each locus there is not a single 
allele that we are counting, but instead we are counting the
number of copies of each haplotype $k\in\{1, 2, \dots, K\}$ at any
specified locus.  Indeed, when $K=2$ we can cast the problem to
SNP case by identifying one of the two haplotypes as an alternative
allele. However when $K>2$ the genetic data is no longer a scalar
but a $K$ dimensional vector $h\in \\{0,1,2\\}^{K\times 1}$ such that
$\sum _{k=1}^K h_k = 2$.  Moreover, this implies that at each locus $j$
there are $K$ effect sizes, that can be expressed as the column
vector $\boldsymbol{\beta}_j\in \mathbb{R}^{K\times 1}$.  As we can
see the haplotype model is more complex as the fixed effects are
in a $K$ dimensional space as opposed to a 1 dimensional space.

An important question when working under the haplotype model is what
type of random genetic effects do we want to account for.  If we
only care about the polygenic SNP effects, then we should make use of
GRMs $\mathbf{A}$ or $\mathbf{A}_\text{EHC}$.  Another choice would
be to account for poly-haplotype effects, that the LMM mapping genotype
to phenotype at locus $j$ for sample $i$ becomes,

$$
\begin{equation}
Y_i = h_{ij}^T\\,\boldsymbol{\beta}_j + \mathbf{W}_1 U_1 + \mathbf{W}_2 U_2 + 
\dots + \mathbf{W}_K U_K + \epsilon_i.
\end{equation}
$$

Here, the difference between the polygenic SNP and haplotype effects
are explicit.  Instead of a single design matrix $\mathbf{Z}$ accounting
for polygenic effects there are now $K$ matrices.  Consequently, the
haplotype GRM is

$$
\mathbf{A}_\text{EHC} = \sum_{k=1}^K \mathbf{W}_k\mathbf{W}^T_k
$$

the sum of the similarity matrices of each haplotype.

### Leave one chromosome out (loco) GRM

In a genome wide association study (GWAS) polygenic effects are often
accounted for from all loci except those of the chromosome that fixed
effect sizes are being estimated.  This approach is referred to as
"leave-one-chromosome-out", or loco for short, for which we will refer
the resulting GRM as the loco GRM.  This matrix can be computed easily
from the set of all chromosome GRMs.

Let $M_c$ be the number of markers used for computing the GRM $A_c$.
Then it follows that the loco GRM $L_u$ that is used for effect size
estimation of all loci on chromosome $u$ is,

$$
\begin{align}
L_u &= \sum_{\forall c \neq u} A_c.
\end{align}
$$


## Command line interface <a name="cli"></a>

The `grm` program consists of two subprograms: 

* `grm contig`: the computation of the GRM of a named contig
* `grm loco`: the aggregation of contig GRMs into a 
    "leave-one-chromosome-out" matrix (denoted the "loco" matrix).  

to read documentation on the respective subprogram simply

```
grm [contig | loco] --help
```

The GRM calculation requires tha the genetic information is in the
the bcf family of file formats, i.e. vcf, vcf.gz, or bcf.  While
the loco subprogram requires all chromosome matrices to be in the
`.grm` binary file format defined below.

### Computing a contig GRM

By default, the GRM is computed
using the expected haplotype counts with FORMAT ID = "HD".

```
grm chrm <chrm_id> <path/to/my/snps.bcf>
```
will produce a binary `.mat` file that stores the GRM and relavent meta data.
Other options include

### Computing a loco GRM




## Installation and requirements <a name="install"></a>

The program is only available as source from this repository and requires

* `GNU make` 
* `htslib` https://github.com/samtools/htslib
* `argparse` https://github.com/robert-vogel/argparse
* `clang` or `gcc` C++17 compiler


## The `.grm` file format <a name="grmspec"></a>

The `.grm` file format is a binary data format consisting of meta data
and a payload.


Assume 64-bit machine, little-endian, e.g. ARM and x86-64.
    

### Meta data
    
| offset    | type      | size      | description                               |
| (bytes)   |           | (bytes)   |                                           |
| --------- | --------- | --------- | ----------------------------------------- |
| 0         | uint32_t  | 4         | File signature (0x47524D00 = "GRM\0")     |
| 4         | uint32_t  | 4         | File version, utils::Version              |
| 8         | uint32_t  | 4         | Program version, utils::Version           |
| 12        | varies    | $n_{coords}$      | Genomic coordinates, grm::Coordinates     |
| $n_{coords}$ + 12 | varies | $n_{samps}$  | Sample names, grm::Samples                |
| $n_{samps} + n_{coords} +12$ | | | Total    |


**Genomic coordinates**

| offset    | type      | size      | description                               |
| (bytes)   |           | (bytes)   |                                           |
| --------- | --------- | --------- | ----------------------------------------- |
| 0             | uint32_t  | 4         | $l_{contig}=$ strlen(contig name)         |
| 4 | char[]    | $n_{contig} = l$      | Contig name chars, sizeof(char) = 1 byte  |
| $n_{contig} + 4$ | uint32_t | 4  | Number of positions ($l_{pos}$)           |
| $n_{contig} + 8$ | uint32_t | $n_{pos} = 4 l_{pos}$ | loci positions on contig |
| $n_{coords} =n_{contig} + n_{pos} + 8$  | |          Total                 |


**Sample names**

| offset    | type      | size      | description                               |
| (bytes)   |           | (bytes)   |                                           |
| --------- | --------- | --------- | ----------------------------------------- |

### Payload

Given $N$ samples, the payload is the upper triangular components of the
$N \times N$ GRM stored as an array of $N(N+1)/2$ 32-bit floating point
numbers in row-major order. While the genotype-based GRM will not produce
fractional values, the expected count GRMs will, making `float32` an
acceptable choice for all GRM types.



## Contributing <a name="contributing"></a>

I am using [GoogleTest](https://google.github.io/googletest/) framework
for organizing tests.  If you contribute, please make tests for your
contributions.  To run tests, `build` directory and build the project
make check



## A.I. Acknowledgement <a name="ai"></a>

The problem statement and overall design of the code base was by 
Robert Vogel. He has made use of Claude for review and Claude Code,
the AI assistant by Anthropic, to implement sum features.


## References <a name="refs"></a>

[1] [Yang et al. Nature Genetics 42, 565-569 (2010)](https://www.nature.com/articles/ng.608)
[2] [Davies et al. Nature Genetics 48, 965-969 (2016)](https://www.nature.com/articles/ng.3594)
<!--
[3] [Kang et al. Genetics 178: 1709-1723 (2008)](https://academic.oup.com/genetics/article/178/3/1709/6061473)

[2] [Kang et al. Nature Genetics 42 348-354 (2010)](https://www.nature.com/articles/ng.548)
-->
