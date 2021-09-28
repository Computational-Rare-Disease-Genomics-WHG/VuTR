""" Finds the uORFs within a sequence"""

import pandas as pd
from utr_utils.tools.mane import read_mane_transcript
from utr_utils.tools.mane import find_uorfs_in_transcript

write_location = "../../data/pipeline/UORFS.tsv"


def main():
    # Read in mane

    summary =  # read summary mane

    uorfs_df = pd.DataFrame()
    for key, row in summary.iterrows():
        print(f'Finding uORFS in {row["Ensembl_nuc"]}')
        transcript = read_mane_transcript(row['Ensembl_nuc'])
        uorfs = find_uorfs_in_transcript(transcript)
        append_row = pd.DataFrame(uorfs)
        uorfs_df = pd.concat([uorfs_df, append_row])

    # Write file
    uorfs_df.to_csv(write_location, sep="\t")


if __name__ == '__main__':
    main()
