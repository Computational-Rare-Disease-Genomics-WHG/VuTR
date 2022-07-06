"""Creates a TSV file that creates all possible UTR variants for VEP
TODO : Add indels as well
"""

from itertools import chain, product
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

bases = ['A', 'C', 'G', 'T']
ASSEMBLY = 'GRCh38'
complement_bases = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

def generate_nbases(n: int):
    """
    Generates an array of insertions for a given length
    @param n int : n >= 2
    @returns possible_vars: list of strings 
    """
    possible_vars = [] 
    for i in range(2, n+1):
        possible_vars = possible_vars + [''.join(comb) for comb in product(bases, repeat=i)]
    return possible_vars

def get_reverse_complement(seq):
    """
    Flips the string to it's reverse complement

    @param seq : str 
    @returns reverse_complement : str
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def variant_id(x, seq, start):
    """
    Generates the variant id given the columns

    @param x: Dataframe series with chr, pos and ref attributes 
    @param seq: The 5' UTR sequence 
    @param start: the start coordinate
    @returns variant_id: str the variant id as per gnomad 
    """
    prebase = seq[x['pos']-start-1:x['pos']-start]
    return f'{x["chr"]}-{x["pos"]}-{prebase+x["ref"]}-{prebase}'

def generate_deletions(sequence, chr, start, end, strand, n):
    """
    Generates a dataframe of deletions of size <= n

    @param sequence : 5' UTR sequence
    @param start: int start of the above sequence
    @param end: int end of the above sequence
    @param strand : char which strand the sequence is located in
    @param n: int the size of the deletion
    @returns var : Dataframe of simulated deletions for the sequence 

    TODO: Add the padding for the boundary cases of the regions
    """
    # Get the reverse complement 
    if strand == '-':
        sequence = get_reverse_complement(seq=sequence)
    
    variant_df = []
    for i in range(n+1)[1:]:
        refs = [sequence[x:x + i] for x in range(len(sequence))][0:(len(sequence) if i==1 else -(i-1))]
        nrow = len(refs)
        # Generate the data frame
        var = pd.DataFrame({'chr' : chr, 
                            'pos': list(range(start, end-i+1)),
                            'strand' : strand, 
                            'ref' : refs,
                            'alt' : '-'
                            })
        var['variant_id'] = var.apply(lambda x: variant_id(x, sequence, start), axis=1)
        # Generate the variant ID column
        variant_df.append(var)

    # Concatenate the dataframe 
    return pd.concat(variant_df)


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
    script_path = Path(__file__).parent  # noqa: E501 # pylint: disable=C0301
    mane_version = args.mane_version
    transcript_sequences = pd.read_csv(
        script_path
        / f'../../data/pipeline/MANE/{mane_version}/MANE_transcripts_v{mane_version}.tsv',  # noqa: E501 # pylint: disable=C0301
        sep='\t',
    )  # noqa: E501 # pylint: disable=C0301
    features = pd.read_csv(
        script_path
        / f'../../data/pipeline/MANE/{mane_version}/MANE.GRCh38.v{mane_version}.ensembl_genomic.tsv',  # noqa: E501 # pylint: disable=C0301
        sep='\t',
    )
    features = features[features['type'] == 'five_prime_UTR']
    features['width'] = features['end'] - features['start'] + 1

    # Write output
    features = features.sort_values(by=['gene_id', 'exon_number'])
    chroms = list(range(1, 23)) + ['X', 'Y']

    if args.only_chr_22:
        # filter utr_file to chr_22
        chroms = ['22']

    if not args.include_indels : 
        insertion_array = generate_nbases(args.indel_size)

    formated_chroms = ['chr' + str(i) for i in chroms]
    
    for chrom in formated_chroms:
        vprint(f'Starting generating mutations for {chrom}')
        chrom_possible_df = pd.DataFrame()
        long_df_list = []
        for gene in features[features['seqid'] == chrom]['gene_id'].unique():

            # find all of utr features  within that region
            feats = features[features['gene_id'] == gene]
            chrom = feats['seqid'].values[0]
            strand = feats['strand'].values[0]
            utr_length = sum(feats['width'])

            # find the sequence
            seqs = list(
                transcript_sequences[transcript_sequences['gene'] == gene]
                .seq.values[0][0:utr_length]
                .upper()
            )  # pylint: disable=C0301
            # Get the list of positions
            pos = list(
                chain(
                    *list(
                        feats.apply(
                            lambda x: list(range(x['start'], x['end'] + 1)), axis=1
                        )
                    )
                )
            )  # noqa: E501 # pylint: disable=C0301

            if strand == '-':
                # sort based on strand
                pos = sorted(pos, reverse=True)

                # get complement of the sequence
                seqs = [complement_bases[nt] for nt in seqs]

            # Simulate a dataframe of all variants
            long_df = pd.DataFrame(
                {
                    'chrom': [chrom[3:]] * (utr_length * 4),
                    'start': np.repeat(pos, 4),
                    'end': np.repeat(pos, 4),
                    'ref': np.repeat(seqs, 4),
                    'alt': bases * utr_length,
                    'strand': [strand] * (utr_length * 4),
                }
            )


            # Generate all insertions and add them up
            if args.include_indels:
                long_df = pd.concat([long_df, 
                    pd.DataFrame(
                        {
                            'chrom' : [chrom[3:]] * utr_length * len(insertion_array),
                            'start' : np.repeat(pos, len(insertion_array)),
                            'end' : np.repeat(pos, len(insertion_array)),
                            'ref' : np.repeat(seqs, len(insertion_array)),
                            'alt' : insertion_array * utr_length,
                            'strand' : [strand] * (utr_length*len(insertion_array))
                        }
                    )
                ])

                # Generate all deletions 
                

            break
            # Remove all rows with ref same as alt
            long_df = long_df[long_df['ref'] != long_df['alt']]

            # Format to VEP input
            long_df['allele'] = long_df['ref'] + '/' + long_df['alt']
            long_df = long_df.loc[:, ['chrom', 'start', 'end', 'allele', 'strand']]
            long_df_list.append(long_df)
        break
        chrom_possible_df = pd.concat(long_df_list, ignore_index=True)

        vprint(f'Finish generating mutations for {chrom}')
        vprint(f'Writing to VEP file', args.verbose)
        chrom_possible_df = chrom_possible_df.sort_values(by='start').drop_duplicates()
        chrom_possible_df.to_csv(
            script_path
            / f'../../data/pipeline/vep_data/input/UTR_variants_all_possible_{ASSEMBLY}_{mane_version}_{chrom}.txt',  # noqa: E501 # pylint: disable=C0301
            sep='\t',
            header=None,
            index=False,
        )
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates a TSV file that creates all possible UTR variants'  # noqa: E501 # pylint: disable=C0301
    )
    parser.add_argument(
        '--mane_version',
        default='1.0',
        help='Which mane_version to use?',
    )
    parser.add_argument(
        '--indel_size', 
        default=3, 
        help='How large should the simulated indels be?'
    )
    parser.add_argument(
        '--include_indels',
        action='store_true',
        default=False,
        help='Whether to exclude indels of (2bps), defaults to false',
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
