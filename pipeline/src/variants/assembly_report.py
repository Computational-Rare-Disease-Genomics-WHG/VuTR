
import argparse
import pandas


def main(args):
    """
    Entrypoint
    """
    b = pandas.read_csv(
        args.input,
        sep="\t",
        skiprows=61)
    b = b[['# Sequence-Name', 'RefSeq-Accn', 'UCSC-style-name', 'Assembly-Unit']]
    b = b[b['Assembly-Unit'] == "Primary Assembly"]
    b = b[b['# Sequence-Name'].isin([str(i) for i in list(range(1, 23))+["X", "Y"]])]
    b.to_csv(args.output, index=False, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=""  # noqa: E501 # pylint: disable=C0301
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input File",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file",
    )
    main(args=parser.parse_args())
