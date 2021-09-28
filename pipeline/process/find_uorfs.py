""" Finds the uORFs within a sequence
TODO: Run for each transcript
"""

import pandas as pd
from utr_utils.tools.mane import read_mane_transcript


def main():
    # Read in mane

    summary =
    mane = read_mane_transcript()

    print(mane)

    for key, row in mane.iterrows():

        find_uorfs_in_transcript()
        pass
    pass


if __name__ == '__main__':
    main()
