# Process MANE

Set of scripts to run to process MANE files

## TODO

- [ ] Get rid of hardcoded path
- [ ] Fix working directories
- [ ] Generate regions sort of like (`find_utr_regions.R`) files for 15 nt introns.
- [ ] (Optional) Rewrite in python to be consistent with the rest

## Running this set of scripts

```bash
# Writes the mane gff file to tsv
Rscript convert_mane_features_to_tsv.R

# Creates a file demarcating UTR genomic locations
Rscript find_utr_regions.R

# Parses the transcript files and annotates with CDS/UTR features
Rscript transcripts.R

# Creates a table mapping genomic locations with transcript locations
Rscript create_gpos_lookup.R

# Finds the uORFS and relatevant fetures per transcript id
Rscript find_uorfs.R

# Finds the genomic locations of the uORFS
Rscript find_genomic_locations.R
```
