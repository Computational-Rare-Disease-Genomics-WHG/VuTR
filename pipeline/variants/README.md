# Variants


- Generates a list of possible variants (SNVs inside UTRs for now)
- Subsets gnomAD and 

## Future log

- [ ] merge_possible_variants.py
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