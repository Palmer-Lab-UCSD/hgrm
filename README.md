# 🏗️ Being built 🏗️

# Compute the genetic relationship matrix


The genetic relationship matrix (GRM) describes the genetic relationship 
between pairs of samples.  Its computation is dependent on the random effects
defined by the linear mixed model mapping genetic features to phenotype.  For
example, suppose that we are interested in accounting for polygenic SNP effects
using measured genotypes.  Let the number of samples be $N$, the number of
loci genotyped $M+1$, $Y\in\mathbb{R}^{N\times 1}$

$$
\begin{align}
Y &= x_j\beta + \mathbf{Z}_j U_j + \epsilon\\
U_j &\sim \mathcal{N}\left(0, \sigma_g^2 \mathbf{I}_{M\times M}\right)\\
\epsilon &\sim \mathcal{N}\left(0, \sigma_e^2 \mathbf{I}_{N\times N}\right)
\end{align}
$$

genotype vector at locus $j$, $\mathbf{Z}_j$ is an $N\times M$ genotype matrix
consisting of a set loci that do not include locus $j$, 
$U_j\sim\mathcal{N}\left(0,\sigma^2_g \mathbf{I}_{M\times M}\right)$
independent genetic random effects, and 
$\epsilon\sim\mathcal{N}\left(0,\sigma^2_e\mathbf{I}_{N\times N}\right)$
independent environmental random effects.  From which it follows that

$$
\text{cov}(Y) = \sigma^2_g\mathbf{Z}\mathbf{Z}^T
+ \sigma^2_e\mathbf{I}_{N\times N}
$$

where the genetic component of the phenotype covariance tells us the how
to compute the GRM, i.e. $\mathbf{Z}\mathbf{Z}^T$.

This program provides the tools to compute the GRM genotypes (alt allele
count 0,1,2), the expected alt allele count under the imputation model,
the expected haplotype count, and a combination of the expected alt allele
count with the expected haplotype count.


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


## Installation and availability

The program is only available as source from this repository and requires

* `GNU make` 
* `htslib`
* `clang` or `gcc` C++17 compiler



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

## References

[1] [Kang et al. Genetics 178: 1709-1723 (2008)](https://academic.oup.com/genetics/article/178/3/1709/6061473)

[2] [Kang et al. Nature Genetics 42 348-354 (2010)](https://www.nature.com/articles/ng.548)

[3] [Yang et al. Nature Genetics 42, 565-569 (2010)](https://www.nature.com/articles/ng.608)

