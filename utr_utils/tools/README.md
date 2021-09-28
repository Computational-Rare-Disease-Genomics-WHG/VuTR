# Tools 


These are the various different tools that the UTR Visualization app will be interfacing with. 

## gnomAD 

We can access gnomAD through the following query through the GraphQL API.  

Can play with it [here.](https://gnomad.broadinstitute.org/api)

```js
{
  transcript(transcript_id: $transcript_id, reference_genome: GRCh38) {
    transcript_id
    transcript_version
    strand
    start
    stop
    gene {
      gene_id
      name
      omim_id
    }
    clinvar_variants {
      transcript_id
      ref
      pos
      alt
      in_gnomad
      clinvar_variation_id
      gold_stars
      review_status
      hgvsc
      clinical_significance
    }
    variants(dataset: gnomad_r3) {
      ref
      pos
      alt
      genome {
        af
        an
        ac
      }
      transcript_consequence {
        is_mane_select
        sift_prediction
        polyphen_prediction
        is_mane_select_version
      }
      variantId
    }
  }
}

```


## UTR Annotator 

We do this through a bash script that is run from a python interface

```
```


## sORFs.org 

Through the biomaRt 


## ClinGen (https://dosage.clinicalgenome.org/)

Through Entrez Utils?