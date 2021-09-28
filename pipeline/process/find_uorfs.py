""" Finds the uORFs within a sequence
TODO: Run for each transcript
"""

import pandas as pd
from utr_utils import find_uorfs_in_transcript


def main():
    # Read in mane
    mane = pd.read_csv()
    for key, row in mane.iterrows():

        find_uorfs_in_transcript()
        pass
    pass


if __name__ == '__main__':
    main()
