"""
generate.py

A command line tool to generate variants from a reference genome.
Generate all possible variants output as a vcf file

Usage:
    python generate.py \
        --reference_fasta_file /path/to/reference.fa \
        --assembly_report /path/to/assembly_report.txt \
        --bed_file /path/to/bed_file.bed \
        --output_file /path/to/output.vcf \
        --indel_size 3 \
        --num_processes 4 \
        --only_chr_22
"""

import argparse
import itertools
import pandas

from pathlib import Path
from multiprocessing import Pool
from tqdm import tqdm
from pyfaidx import Fasta


def generate_snvs(
    sequence: str,
    seq_start_pos: int,
    chrom: str
):
    """
    Generate SNVs given a sequence
    
    @param sequence : The sequence to generate variants from
    @param seq_start_pos : The start position of the sequence
    @param chrom : The chromosome the sequence is on
    @returns variants : A list of tuples
        containing the position, reference, alternate, and variant id
    """
    variants = []
    for i, _s in enumerate(sequence):
        for nt in ["A", "C", "G", "T"]:
            if nt != sequence[i]:
                var_id = f"{chrom}-{seq_start_pos+i-3}-{sequence[i]}-{nt}"
                variants.append((seq_start_pos+i-3, sequence[i], nt, var_id))
    return variants


def generate_insertions(sequence: str, seq_start_pos: int, chrom: str, n: int):
    """
    Generate insertions given the sequence

    @param sequence : The sequence to generate variants from
    @param seq_start_pos : The start position of the sequence
    @param chrom : The chromosome the sequence is on
    @param n : The maximum length of the insertion

    @returns variants : A list of tuples containing the position,
        reference, alternate, and variant id
    """
    variants = []
    for i in range(len(sequence) + 1):
        for ins_len in range(1, n + 1):
            for ins_nts in itertools.product(["A", "C", "G", "T"], repeat=ins_len):
                ins_seq = "".join(ins_nts)
                if ins_seq:
                    if i >= 1:
                        var_id = f"{chrom}-{seq_start_pos + i + 2}-{sequence[i-1]}-{sequence[i-1]+ins_seq}"
                        variants.append((seq_start_pos + i + 2, sequence[i-1], sequence[i-1]+ins_seq, var_id))
    return variants


def generate_deletions(sequence: str, seq_start_pos: int, chrom: str, n: int):
    """
    Generate deletions given the sequence

    @param sequence : The sequence to generate variants from
    @param seq_start_pos : The start position of the sequence
    @param chrom : The chromosome the sequence is on
    @param n : The maximum length of the deletion

    @returns variants : A list of tuples containing the position
        reference, alternate, and variant id
    """
    variants = []
    for i in range(len(sequence)):
        for del_len in range(1, n + 1):
            if i + del_len <= len(sequence):
                del_seq = sequence[i : i + del_len]
                pre_base = sequence[i - 1]
                var_id = f"{chrom}-{seq_start_pos+i - 4}-{pre_base+del_seq}-{pre_base}"
                variants.append((seq_start_pos + i - 4,
                                 pre_base+del_seq, pre_base, var_id))
    return variants


def generate_vcf_lines(chrom, variants):
    """create a set of variants to avoid duplicates"""
    variant_set = set()
    vcf_str = []
    for pos, ref, alt, var_id in variants:
        if ref == alt:
            continue
        if (pos, ref, alt, var_id) in variant_set:
            continue
        # Add to variant set
        variant_set.add((pos, ref, alt, var_id))
        pos_str = str(pos)
        var_line = f"{chrom}\t{pos_str}\t{var_id}\t{ref}\t{alt}\t{'.'}\t{'.'}\t{'.'}\n"
        vcf_str.append(var_line)
    return vcf_str


def write_vcf_header(output_path):
    """
    Write the vcf header to a file

    @param output_path : The path to write the vcf header to
    @returns None
    """
    header = ["##fileformat=VCFv4.3\n"]
    header.append("##reference=GRCh38\n")
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    with open(output_path, "w", encoding="utf-8") as file:
        file.writelines(header)


def write_to_vcf(vcf_lines, output_path):
    """
    Write the vcf lines to a vcf file

    @param vcf_lines : A list of vcf lines
    @param output_path : The path to write the vcf file to
    @returns None
    """
    with open(output_path, "a", encoding="utf-8") as file:

        # Write VCF lines
        file.writelines(vcf_lines)


def gen_variants(
    sequence: str, 
    seq_start_pos: int, n: int, chrom: str, num_processes=1
):
    """
    Generate variants parallely

    @param sequence : The sequence to generate variants from
    @param seq_start_pos : The start position of the sequence
    @param n : The maximum length of the insertion/deletion
    @param chrom : The chromosome the sequence is on
    @param num_processes : The number of processes to use
    @returns variants : A list of tuples containing the position,

    """
    with Pool(num_processes) as pool:
        snvs = pool.apply_async(generate_snvs, [sequence, seq_start_pos, chrom])
        insertions = pool.apply_async(
            generate_insertions, [sequence, seq_start_pos, chrom, n]
        )
        deletions = pool.apply_async(
            generate_deletions, [sequence, seq_start_pos, chrom, n]
        )
        variants = snvs.get() + insertions.get() + deletions.get()

    return variants


def main(args):
    """
    Main function, command line arguments are passed as args
    """
    reference = Fasta(args.reference_fasta_file)

    chrom_set = [str(i) for i in list(range(1, 23)) + ["X", "Y"]]
    if args.only_chr_22:
        chrom_set = "22"

    assembly_report = pandas.read_csv(args.assembly_report, sep="\t")

    with open(args.bed_file, encoding="utf-8") as bedfile:
        # Filter bed file to chromosome
        for line in tqdm(bedfile):
            
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            write_vcf_header(args.output_file)
            if chrom in chrom_set:                
                # Get sequence for region
                acc_id = assembly_report[
                    assembly_report['# Sequence-Name'] == chrom]["RefSeq-Accn"].values[0]
                exon_seq = reference[acc_id][start - (args.indel_size+1) : end].seq.upper()

                # Generate variants from the sequence
                generated_vars = gen_variants(
                    exon_seq, start, args.indel_size, chrom, num_processes=args.threads
                )
                # Generate VCF lines
                vcf_lines = generate_vcf_lines(chrom, generated_vars)

                # Write to VCF file
                write_to_vcf(vcf_lines, args.output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Creates a TSV file that creates all possible UTR variants"  # noqa: E501 # pylint: disable=C0301
    )
    parser.add_argument(
        "--bed_file",
        type=str,
        required=True,
        help="Which mane_version to use?",
    )
    parser.add_argument(
        "--reference_fasta_file",
        type=str,
        required=True,
        help="Which mane_version to use?",
    )

    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Which mane_version to use?",
    )

    parser.add_argument(
        "--assembly_report",
        type=str,
        required=True,
        help="Which mane_version to use?",
    )

    parser.add_argument(
        "--indel_size",
        type=int,
        default=3,
        help="How large should the simulated indels be?",
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="How many threads to use?"
    )
    parser.add_argument(
        "--include_indels",
        action="store_true",
        default=False,
        help="Whether to exclude indels of (3bps), defaults to false",
    )
    parser.add_argument(
        "--only_chr_22",
        action="store_true",
        help="Verbose outputs",
    )
    main(args=parser.parse_args())
