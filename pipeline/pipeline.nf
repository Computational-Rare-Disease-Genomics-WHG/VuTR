/*
* Converts to TSV
*/
process convertToTsv{
    input: 
    file input_file from '{$params.baseDir}/{$params.data_directory}/MANE/{$params.MANE.mane_version}/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz'

    output:
    file mane_gff.tsv to '{$params.baseDir}/{$params.data_directory}/MANE/{$params.MANE.mane_version}/'

    script: 
    """
    Rscript src/process_mane/convert_mane_features_to_tsv.R \
         --mane_gff {$input_file} --output_file ${}
    """
}


process tableTranscripts{
    input: 

    output: 

    script: 
    """
    Rscript src/process_mane/transcripts.R --mane-file ${} --rna-file {} --output-rna-file ${} --output-rna-feature-file ${}
    """
}

process createVariantDatabase{
    input: 

    output: 

    script: 
    """
    """
}

process runVEP{
    input:
    file "input.txt"

    output: 
    file "output.txt"

    containerOptions {

    }


    switch (params.profile){
        case "local": 
            script: 
            """
            vep \
            --assembly GRCh38\
            --force_overwrite\
            --species homo_sapiens\
            --cache \
            --fasta ${VEP_CACHE}/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
            --mane \
            --mane_select\
            --canonical\
            --offline\
            --tab \
            --fork 10\
            --dir_cache ${VEP_CACHE} \
            --plugin UTRannotator\
            -i ${INPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr${i}.txt\
            -o ${OUTPATH}/UTR_variants_vep_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr${i}.txt
            """


        case "cluster":
            script: 
            """
            vep \
            --assembly GRCh38\
            --force_overwrite\
            --species homo_sapiens\
            --cache \
            --fasta ${VEP_CACHE}/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
            --mane \
            --mane_select\
            --canonical\
            --offline\
            --tab \
            --fork 10\
            --dir_cache ${VEP_CACHE} \
            --plugin UTRannotator\
            -i ${INPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr${i}.txt\
            -o ${OUTPATH}/UTR_variants_vep_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr${i}.txt
            """
    }


    
}


process createFeatureDatabase{
    input: 

}





workflow {
    processA()
    processB()
    processC()
    
    output:
    file "final_output.txt" into params.output_dir
}







// /*w
// Creates a genomic lookup table 
// */
// process createLookupTable {
//     input:
//     file input_file from convertToTsv
//     file rna_file 

//     output:
//     file dx
//     file file

//     script:
//     """
//     Rscript src/process_mane/create_gpos_lookup.R
//     """
// }


// // Create smORFS
// process manageSmORFs{
//     input:
//     file fxc from
//     file cx from

//     output:

//     script:
//     """
//     Rscript  \
//     --mane-file  {}\
//     --smorf {} \
//     --g2tcoord {} \
//     --output-file {}\
    
//     """
// }



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

// // process variantProcess{

// //     // Replace run_vep_parser_on_all to this file
// //     script:
// //     """
// //     python3 generate_all_possible_variants.py --include_indels
// //     run_vep_parser_all_chrs.sh
// //     """

// // }

// // process ingestFeaturesDatabase{
// //     // input : 

// //     script : 
// //     """
// //     python3 ingest_table.py --db_name features.db --all --verbose --overwrite
// //     python3 variant_store.py --db_name variant_store.db --variant_file 
// //     """


// // }
