### SarahSeq

1. Install the R package with `devtools::install_github('Sarah145/SarahSeq')`
2. Load the library with `library(SarahSeq)`
3. Run the shiny app with `SarahSeq()`
4. When app loads, upload your snp.vcf.gz file from Dante Labs (must be gzipped, or else you'll get max upload file size exceeded error)
5. Play!

Disclaimer: This whole process is very slow, especially step 1 because the package includes uk biobank data and step 4 because your vcf contains ~4 million variants. I'm not sure if there is a way to optimise any of the steps because the data is just so large but I'm open to suggestions!



