Please read [here](https://weinstockj.github.io/htslibr) for documentation of the htslibr R package. To install
this package, you can use 
```git clone --depth=50 --branch=master https://github.com/weinstockj/htslibr.git weinstockj/htslibr
cd weinstockj/htslibr/
git submodule update --init --recursive
cd htslibr/
R CMD build .
```

This package was built on a Ubuntu 16.04 system, and has not been tested for portability to 
OS X or Windows. Not likely to work anywhere where htslib is difficult to compile. 

## References

This R package is bundling the C htslib library. 
[htslib](https://github.com/samtools/htslib)

Previous efforts to wrap htslib and build functionality include the following packages. These packages
are more mature than htslibr. 
[Rhtslib](http://bioconductor.org/packages/release/bioc/html/Rhtslib.html)
[Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)

