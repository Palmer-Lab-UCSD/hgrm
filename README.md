# 🏗️ Being built 🏗️

# Compute the genetic relationship matrix using expected haplotype counts


The genetic relationship matrix (GRM) describes the genetic relationship between
pairs of samples.  WRITE MORE




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

