
## Background 
The 5’ UTR of a given gene includes cis-elements that regulate mRNA translation such as uORFs. Variants in 5’ UTRs can modify or entirely disrupt translation through various different mechanisms (such as new uAUG that can create uORFs or ORF perturbing variants). Previous work has characterized the role of these in various different loss-of-function mechanisms [Whiffin et al. 2020, Wright et al. 2021] and association in genetic disease. Much of this work has been built using UTR Annotator (https://github.com/ImperialCardioGenetics/UTRannotator) , which is a VEP plugin built in perl that is able to annotate the consequence of high-impact 5’ UTR variants. 

This project improves accessibility of the UTR Annotator to different researchers through a freely accessible, and user-friendly web app who (i) may not have expertise working with VEP / command-line utility tools, (ii) want a quick method of browsing through a certain gene 5’ UTR or (iii) want to identify allele frequencies in the general population (iv) search if their gene / variant of interest has been classified and reviewed in OMIM / sORFs.org / ClinVar as a potential deleterious gene-of-interest. 


## Components 

### Pipeline
A set of scripts are available at `pipeline`, these are preprosessing scripts that will identify all 5' UTR Variants that are deleterious (in ClinVar) and population variants (in gnomAD) and then run UTR Annotator on them. The outcome will be a Sqlite3 database that will store variants with the consequence on the MANE transcript. 

### Web server 

Source code for the web application that will visualize the 5' UTR Architchure. Currently this is hosted on a Hertzner Server using Docker in [http://49.12.238.72:8080/](http://49.12.238.72:8080/).
