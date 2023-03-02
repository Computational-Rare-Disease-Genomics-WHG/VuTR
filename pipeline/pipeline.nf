
process maneProcess{
    input :  
    output :


    script : 
    """
    Rscript convert_mane_features_to_tsv.R

    # Creates a file demarcating UTR genomic locations
    Rscript find_utr_regions.R

    # Parses the transcript files and annotates with CDS/UTR features
    Rscript transcripts.R

    # Creates a table mapping genomic locations with transcript locations
    Rscript create_gpos_lookup.R

    # Finds the uORFS and relatevant fetures per transcript id
    Rscript find_orfs.R

    # Finds the genomic locations of the uORFS
    Rscript find_genomic_locations.R

    # Find smorfs
    Rscript smorfs.R
    """
}

process variantProcess{

    // Replace run_vep_parser_on_all to this file
    script:
    """
    python3 generate_all_possible_variants.py --include_indels
    run_vep_parser_all_chrs.sh
    """

}

process ingestFeaturesDatabase{
    // input : 

    script : 
    """
    python3 ingest_table.py --db_name features.db --all --verbose --overwrite
    python3 variant_store.py --db_name variant_store.db --variant_file 
    """


}

