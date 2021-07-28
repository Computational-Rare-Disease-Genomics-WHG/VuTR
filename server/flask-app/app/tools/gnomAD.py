import json
import requests

# From playing around in 
# https://gnomad.broadinstitute.org/api?operationName=getVariant&query=%7B%0A%20%20transcript(transcript_id%3A%20%22ENST00000340208%22%2C%20reference_genome%3A%20GRCh38)%20%7B%0A%20%20%20%20transcript_id%0A%20%20%20%20transcript_version%0A%20%20%20%20strand%0A%20%20%20%20start%0A%20%20%20%20stop%0A%20%20%20%20gene%20%7B%0A%20%20%20%20%20%20gene_id%0A%20%20%20%20%20%20name%0A%20%20%20%20%20%20omim_id%0A%20%20%20%20%7D%0A%20%20%20%20clinvar_variants%20%7B%0A%20%20%20%20%20%20transcript_id%0A%20%20%20%20%20%20ref%0A%20%20%20%20%20%20pos%0A%20%20%20%20%20%20alt%0A%20%20%20%20%20%20in_gnomad%0A%20%20%20%20%20%20clinvar_variation_id%0A%20%20%20%20%20%20gold_stars%0A%20%20%20%20%20%20review_status%0A%20%20%20%20%20%20hgvsc%0A%20%20%20%20%20%20clinical_significance%0A%20%20%20%20%7D%0A%20%20%20%20variants(dataset%3A%20gnomad_r3)%20%7B%0A%20%20%20%20%20%20ref%0A%20%20%20%20%20%20pos%0A%20%20%20%20%20%20alt%0A%20%20%20%20%20%20genome%20%7B%0A%20%20%20%20%20%20%20%20af%0A%20%20%20%20%20%20%20%20an%0A%20%20%20%20%20%20%20%20ac%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20transcript_consequence%20%7B%0A%20%20%20%20%20%20%20%20is_mane_select%0A%20%20%20%20%20%20%20%20sift_prediction%0A%20%20%20%20%20%20%20%20polyphen_prediction%0A%20%20%20%20%20%20%20%20is_mane_select_version%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20variantId%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D%0A




QUERY = """
query get_data($transcript_id: String!){
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
"""

response = requests.post(
   "https://gnomad.broadinstitute.org/api",
   data=json.dumps({
      "query": QUERY,
      "variables": {"transcript_id": "ENST00000340208"},
   }),
   headers={
      "Content-Type": "application/json",
   },
).json()
print(response)
