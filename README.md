
# VuTR : Visualising the impact of 5' UTR Variants 


![CI Build](https://github.com/Computational-Rare-Disease-Genomics-WHG/VuTR/actions/workflows/tests.yaml/badge.svg)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/Computational-Rare-Disease-Genomics-WHG/VuTR/issues)


## Background 

VuTR is a web application that provides an easy and accessible way to annotate high-impact 5' UTR variants. These variants can modify or disrupt mRNA translation through various mechanisms, such as creating upstream open reading frames (uORFs) or perturbing the ORF. This can have implications for genetic disease and loss-of-function mechanisms.

Previous work has been done using UTR Annotator, a VEP plugin built in Perl that annotates the consequence of high-impact 5' UTR variants. However, this tool can be challenging for researchers who do not have experience with VEP or command-line utility tools. VuTR aims to make this process more user-friendly and accessible.

## Features

VuTR provides the following features:

- An intuitive web interface for easy variant annotation
- The ability to browse through a specific gene's 5' UTR
- Identification of allele frequencies in the general population
- Searching to see if a gene or variant of interest has been classified and reviewed in OMIM, sORFs.org, or ClinVar as a potential deleterious gene-of-interest.

## How to use VuTR

To use VuTR, simply navigate to the website at [vutr.rarediseasegenomics.org](vutr.rarediseasegenomics.org). Once on the homepage, enter your gene as a ENSG or HGNC gene symbol in the search bar. 


## Contributing

Contributions to VuTR are always welcome! Please see the CONTRIBUTING.md file for more information. 


## Installation and Reproducibility

VuTR features two key components, a Nextflow pipeline to generate the databases and a Flask Web server that is used to view vutr.rarediseasegenomics.org. Visit `pipeline/` or `server/` to see additional information on how to install or reproduce VuTR.


## License

VuTR is released under the GPLv2 licencce. Please see the LICENSE.md file for more information.

## Contact

If you have any questions, concerns, or feedback, please reach out to [us](nwhfffin@well.ox.ac.uk). We are always happy to hear from our users and strive to provide the best possible user experience. If you find a bug within the program or have a new feature request, you can alternatively open an Issue on this Github repository. 





