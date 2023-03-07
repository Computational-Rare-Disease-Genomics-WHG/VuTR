// Import the YAML library
import org.yaml.snakeyaml.Yaml

// Load the configuration from the YAML file
def config = new Yaml().load(new File('../config/config.yml').text)


/*
* Converts to TSV
*/
process convertToTsv{
    input: 
    file input_file from config.output

    output: 
    file mane_gff.tsv

    script: 
    """
    Rscript src/process_mane/convert_mane_features_to_tsv.R --mane_gff ${} --output_file ${}
    Rscript src/process_mane/transcripts.R --mane-file ${} --rna-file {} --output-rna-file ${} --output-rna-feature-file ${}
    """
}

/*w
Creates a genomic lookup table 
*/
process createLookupTable {
    input:
    file input_file from convertToTsv
    file rna_file 

    output:
    file dx
    file file

    script:
    """
    Rscript src/process_mane/create_gpos_lookup.R
    """
}


// Create smORFS
process manageSmORFs{
    input:
    file fxc from
    file cx from

    output:

    script:
    """
    Rscript  \
    --mane-file  {}\
    --smorf {} \
    --g2tcoord {} \
    --output-file {}\
    
    """
}


















// process maneProcess{
//     input :  
//     output :


//     script : 
//     """
//     Rscript convert_mane_features_to_tsv.R

//     # Creates a file demarcating UTR genomic locations
//     Rscript find_utr_regions.R

//     # Parses the transcript files and annotates with CDS/UTR features
//     Rscript transcripts.R

//     # Creates a table mapping genomic locations with transcript locations
//     Rscript create_gpos_lookup.R

//     # Finds the uORFS and relatevant fetures per transcript id
//     Rscript find_orfs.R

//     # Finds the genomic locations of the uORFS
//     Rscript find_genomic_locations.R

//     # Find smorfs
//     Rscript smorfs.R
//     #;/bin/bash

// # Get bedtools working
// # bedtools getfasta -fi GRCh38_latest_genomic.fna -bed my_intervals.bed -fo my_sequences.fa

//     """
// }

// process variantProcess{

//     // Replace run_vep_parser_on_all to this file
//     script:
//     """
//     python3 generate_all_possible_variants.py --include_indels
//     run_vep_parser_all_chrs.sh
//     """

// }

// process ingestFeaturesDatabase{
//     // input : 

//     script : 
//     """
//     python3 ingest_table.py --db_name features.db --all --verbose --overwrite
//     python3 variant_store.py --db_name variant_store.db --variant_file 
//     """


// }
