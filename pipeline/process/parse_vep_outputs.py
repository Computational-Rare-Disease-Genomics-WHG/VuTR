"""
Parses the outputs from VEP running 
UTR annotator and saves it as a .tsv
TODO : Add gnomad when getting access to it
"""

import pandas as pd
import numpy as np


def wide_to_long(df):
    """
    Make from wide to long. 
    @param df with 5'utr consequence concatentated by &.
    @returns long_df What we want is a 
            row per variant / transcript / five_prime_UTR_variant_annotation. 
    """
    long_df = pd.DataFrame(columns=df.columns.values.tolist())

    for index, row in df.iterrows():

        binding_df = pd.DataFrame()
        # If single just row bind it to long_df
        if "&" not in row["five_prime_UTR_variant_consequence"]:
            binding_df = pd.DataFrame(row.to_frame().T)
        else:
            # Replicate the rows per ORF consequences
            nrows = row["five_prime_UTR_variant_consequence"].count("&")
            consequences_split = row["five_prime_UTR_variant_consequence"].split("&")
            annotation_split = row["five_prime_UTR_variant_annotation"].split("&")
            binding_df = pd.concat([row.to_frame().T.drop(
                columns=["five_prime_UTR_variant_consequence",
                         "five_prime_UTR_variant_annotation"])]*(nrows+1),
                ignore_index=True)
            binding_df["five_prime_UTR_variant_consequence"] = consequences_split
            binding_df["five_prime_UTR_variant_annotation"] = annotation_split
        long_df = pd.concat([long_df, binding_df], axis=0, ignore_index=True)

    return long_df


# TODO : Update directories to cmd line args
clinvar_vep_path = '../../data/pipeline/vep_data/output/clinvar_utr_filtered.out.txt'
mane_summary_path = '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.summary.txt.gz'
write_path = "../../data/pipeline/CLINVAR/clinvar_processed.txt"
# Read clinvar file and the mane summary file
clinvar_df = pd.read_csv(clinvar_vep_path,
                         sep='\t',
                         skiprows=44)

mane_summary_df = pd.read_csv(mane_summary_path,
                              sep='\t')

# Filter to variants that impact uORFs (remove all other variants)
no_utr_consequence = ['-', '', None, np.nan]
clinvar_df = clinvar_df[~clinvar_df['five_prime_UTR_variant_annotation'].isin(
    no_utr_consequence)]

# Remove version number
mane_summary_df['transcript_id'] = mane_summary_df['Ensembl_nuc'].apply(
    lambda x: x[0:15])

print(mane_summary_df["transcript_id"])

# Filter to consequence on the MANE transcript
clinvar_df = clinvar_df[clinvar_df['Feature'].isin(mane_summary_df['transcript_id'])]


# Convert wide to long data frame
long_clinvar_df = wide_to_long(clinvar_df)

# Write to file
long_clinvar_df.to_csv(write_path,
                       sep="\t",
                       index=False)
