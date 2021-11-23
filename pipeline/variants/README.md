# Variants


- Generates a list of possible variants (SNVs inside UTRs for now, but eventually moving to indels)
- Runs UTR annotator for this list of possible variants
- Subsets gnomAD and Clinvar for UTR variants


## Future log

- [ ] Update directories and I/O for scripts here
- [ ] Generate possible of indels 1-3 nt and then run them through VE!P
- [ ] Expand variant search space to possible variants inside introns & Add VE!P flag for splice variants
- [ ] Append all UTR variants and clinvar annotated variants into gnomAD
- [ ] Convert pandas df in generate_all_possible_variants.py to dask

## Running

```bash 
# Generates all possible variants
python3 generate_all_possible_variants.py --verbose --exclude_indels

# Parse VEP on all shareded possible variants
# Filters out varaints with no impact
bash run_vep_parse.sh

# Merge all possible variants chrom chunks to a single file
python3 merge_possible_variants.py

# Subsets both clinvar and gnomAD
bash run_tabix_clinvar.sh
bash run_tabix_gnomAD.sh

# Extract relevant fields from gnomAD and clinVar
python3 subset_vcfs.py --dataset gnomad
python3 subset_vcfs.py --dataset clinvar
```

## VEP (Docker)

We can run VEP through docker (in environments that can support it). 

But due to compute this may not be possible in instances where we need to use the HPC / BMRC where docker isn't available.

Might need to look into running docker through singularity.

We can interactivately look at the Ensembl docker image using.

```sh
docker pull ensemblorg/ensembl-vep

# To run vep standalone
docker run -t -i ensemblorg/ensembl-vep ./vep

# To view the image interactively running vep
docker run -t -i -v $(realpath ../../data/pipeline/vep_data/)/vep_data/:/opt/vep/.vep/ ensemblorg/ensembl-vep bash
```