# Pipeline 


This is a description of the pipeline to annotate all of the ClinVar and gnomAD UTR variants for the UTR Visualization app. 

To run this pipeline we require docker and the cache. 


```sh 
docker pull ensemblorg/ensembl-vep

# To run vep standalone 
docker run -t -i ensemblorg/ensembl-vep ./vep

# To view the image interactively running vep 
docker run -t -i -v $(pwd)/UTRannotator/:/opt/vep/.vep/Plugins \
 ensemblorg/ensembl-vep \
 bash


# To run vep (with UTR annotator running alongside)
docker run -t -i -v $(pwd)/UTRannotator-vep-plugin:/opt/vep/.vep/Plugins 



docker run -t -i -v /Users/elston.dsouza/Projects/UTR-Visualisation-App/pipeline/UTRannotator-vep-plugin/:/opt/vep/.vep/Plugins ensemblorg/ensembl-vep ./vep -i Projects/UTR-Visualisation-App/pipeline/UTRannotator-vep-plugin/test/test_grch37.vcf --tab --database -plugin UTRannotator -o test.output

```


Once VEP has ru
