"""Creates a TSV file that creates all possible UTR variants for VEP.
TODO : Add indels as well
"""

import argparse
import pandas as pd

from itertools import chain
from pathlib import Path

bases = ["A", "C", "G", "T"]
chroms = list(range(1, 23)) + ["X", "Y"]
mane_version = 0.93
assembly = 'GRCh38'
complement_bases = {
    'A': 'T',
    'G': 'C',
    'C': 'G',
    'T': 'A'
}


def vprint(message: str, verbosity: bool = True):
    """
    Wrapper function to print based on verbosity
    @param message : Message to print to STDOUT
    @param verbosity : The arg.verbose flag.
    @returns None
    """
    if verbosity:
        print(message)


def main(args):
    """
    Run the main script
    @param args : Parsed commandline args
    @returns None
    """
    script_path = Path(__file__).parent

    transcript_sequences = pd.read_csv(script_path /
                                       "../../data/pipeline/MANE/0.93/MANE_transcripts_v0.93.tsv", sep="\t")
    features = pd.read_csv(script_path /
                           '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.tsv',
                           sep='\t',
                           )
    features = features[features['type'] == 'five_prime_UTR']
    features['width'] = features['end'] - features['start']+1

    # Write output

    features = features.sort_values(by=['gene_id', 'exon_number'])
    chrom_possible_df = pd.DataFrame()
    if args.only_chr_22:
        # filter utr_file to chr_22
        chroms = ['22']

    formated_chroms = ["chr"+i for i in chroms]

    for chrom in formated_chroms:
        vprint(f'Starting generating mutations for {chrom}')

        for gene in features[features['seqid'] == chrom]['gene_id'].unique():

            # find all of utr features  within that region
            feats = features[features['gene_id'] == gene]
            chrom = feats['seqid'].values[0]
            strand = feats['strand'].values[0]
            utr_length = sum(feats['width'])

            # find the sequence
            seqs = list(
                transcript_sequences[transcript_sequences['gene'] == gene].seq.values[0][0:utr_length].upper())

            # Get the list of positions
            pos = list(
                chain(*list(feats.apply(lambda x: list(range(x["start"], x["end"]+1)), axis=1))))

            if strand == '-':
                # sort based on strand
                pos = sorted(pos, reverse=True)

                # get complement of the sequence
                seqs = [complement_bases[nt] for nt in seqs]

            # Simulate a dataframe of all variants
            long_df = pd.DataFrame({'chrom': [chrom[3:]]*(utr_length*4),
                                    'start': pos*4,
                                    'end': pos*4,
                                    'ref': seqs*4,
                                    'alt': bases * utr_length,
                                    'strand': [strand] * (utr_length*4)})

            # Remove all rows with ref same as alt
            long_df = long_df[long_df['ref'] != long_df['alt']]

            # Format to VEP input
            long_df['allele'] = long_df['ref'] + "/" + long_df['alt']
            long_df = long_df.loc[:, ['chrom', 'start', 'end', 'allele', 'strand']]
            long_df = long_df.reset_index(drop=True)
            chrom_possible_df = pd.concat([chrom_possible_df, long_df])

        vprint(f'Finish generating mutations for {chrom}')
        vprint(f'Writing to VEP file', args.verbose)
        chrom_possible_df.to_csv(script_path /
                                 f"../../data/pipeline/vep_data/input/UTR_variants_all_possible_{assembly}_{mane_version}_{chrom}.txt",
                                 sep="\t",
                                 header=None,
                                 index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates a TSV file that creates all possible UTR variants'
    )
    parser.add_argument(
        '--exclude_indels',
        action='store_true',
        default=True,
        help='Whether to exclude indels of (2bps)',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose outputs',
    )
    parser.add_argument(
        '--only_chr_22',
        action='store_true',
        help='Verbose outputs',
    )
    main(args=parser.parse_args())
