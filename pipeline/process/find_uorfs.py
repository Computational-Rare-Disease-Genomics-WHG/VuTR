""" Finds the uORFs within a sequence"""

from utr_utils.tools.mane import get_transcript_features

write_location = "../../data/pipeline/UORFS.tsv"


def main():
    # Read in mane
    transcript_id = "ENST00000504921.7"
    transcript_feats = get_transcript_features(transcript_id)
    print(transcript_feats["uORF"])


if __name__ == '__main__':
    main()
