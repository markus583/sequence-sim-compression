import os
import argparse
import sys


def read_and_parse_fasta(file_path: str) -> str:
    """Read a FASTA file and return the sequence as a string.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        str: The parsed sequence.
    """
    with open(file_path, "r") as file:
        lines = file.readlines()
        # skip first (header) line and remove newlines
        sequence = "".join(lines[1:]).replace("\n", "")
    return sequence


def compress_sequence(sequence: str, method="gzip", **kwargs) -> int:
    """Compress the sequence using the specified method and return the length of the compressed sequence.

    Args:
        sequence (str): The sequence to compress.
        method (str, optional): The compression method to use. Can be 'gzip' or '7zip'.
            Defaults to 'gzip'.

    Raises:
        ValueError: If the specified method is not supported.

    Returns:
        int: The length of the compressed sequence.
    """
    # encode the sequence as bytes before compressing
    # this is necessary because gzip only accepts bytes
    encoded_sequence = sequence.encode()
    if method == "gzip":
        import gzip

        return len(gzip.compress(encoded_sequence, **kwargs))
    elif method == "7zip":
        import lzma

        # we use the lzma library for 7zip compression
        # One can also use other compression algorithms within 7zip, e.g., lzma2
        # but we use lzma for efficiency
        return len(lzma.compress(encoded_sequence, **kwargs))
    else:
        raise ValueError("Unsupported compression method")


def calculate_alignment_score(
    file_a: str, file_b: str, method="gzip", **kwargs
) -> float:
    """Calculate the alignment score between two sequences.

    Args:
        file_a (str): Path to the first FASTA file, i.e., first organism.
        file_b (str): Path to the second FASTA file, i.e., second organism.
        method (str, optional): The compression method to use. Can be 'gzip' or '7zip'.
            Defaults to 'gzip'.

    Returns:
        float: The alignment score between the two sequences.
    """
    seq_a = read_and_parse_fasta(file_a)
    seq_b = read_and_parse_fasta(file_b)

    # compress the two sequences separately
    c_a = compress_sequence(seq_a, method)
    c_b = compress_sequence(seq_b, method)
    # concatenate the two sequences and compress them together
    c_ab = compress_sequence(seq_a + seq_b, method)

    # score is the difference between the two compression lengths
    score = (c_a + c_b - c_ab) / max(c_a, c_b)
    # a higher score means a better alignment, so it is inverted
    return 1 - score


def main():
    # if no args are passed, use the example files
    # typically, we would raise an error here
    if len(sys.argv) == 1:
        print("Usage: python calculate_alignment.py <file_a> <file_b> [-m <method>]")
        print("Using example files (cat, cattle) and gzip...")
        file_1 = os.path.join("data", "cat.txt")
        file_2 = os.path.join("data", "cattle.txt")
        score = calculate_alignment_score(file_1, file_2)
    else:
        parser = argparse.ArgumentParser(
            description="Calculate alignment score between two sequences."
        )
        parser.add_argument("file_a", type=str, help="Path to the first FASTA file")
        parser.add_argument("file_b", type=str, help="Path to the second FASTA file")
        parser.add_argument(
            "-m",
            "--method",
            type=str,
            default="gzip",
            choices=["gzip", "7zip"],
            help="Compression method (gzip or 7zip)",
        )
        args = parser.parse_args()

        score = calculate_alignment_score(args.file_a, args.file_b, args.method)
    print(f"Alignment score: {score}")
    return score


if __name__ == "__main__":
    main()
