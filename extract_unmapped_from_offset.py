# extract_unmapped_from_offset.py
# Extract unmapped reads from a BAM file by seeking to the last mapped chunk using the BAI index.
# The BAI index is used to find the maximum virtual offset among all mapped regions.
import argparse
import os
import pysam

def read_uint32(f):
    return int.from_bytes(f.read(4), byteorder="little", signed=False)

def read_uint64(f):
    return int.from_bytes(f.read(8), byteorder="little", signed=False)

def get_last_chunk_end(bai_filename):
    """
    Iterates over the BAI file to find the maximum virtual offset (chunk_end)
    among all mapped regions.
    """
    max_vo = 0
    with open(bai_filename, "rb") as bai:
        # Read magic (4 bytes) and number of references (4 bytes)
        bai.read(4)  # skip magic
        n_ref = read_uint32(bai)
        for _ in range(n_ref):
            n_bins = read_uint32(bai)
            for _ in range(n_bins):
                _bin = read_uint32(bai)  # bin number, not used here
                n_chunks = read_uint32(bai)
                for _ in range(n_chunks):
                    # Each chunk: 8 bytes for chunk_beg, 8 bytes for chunk_end
                    _chunk_beg = read_uint64(bai)
                    chunk_end = read_uint64(bai)
                    if chunk_end > max_vo:
                        max_vo = chunk_end
            # Read number of linear index entries and skip them
            n_intv = read_uint32(bai)
            bai.seek(n_intv * 8, os.SEEK_CUR)
    return max_vo

def extract_unmapped_reads_from_offset(bam_file, bai_file, output_bam):
    # Get the maximum virtual offset from the BAI file
    last_vo = get_last_chunk_end(bai_file)
    print(f"Last mapped virtual offset (from BAI): {last_vo}")

    # Open the BAM file and seek to the computed virtual offset.
    with pysam.AlignmentFile(bam_file, "rb") as inbam:
        # Seek to the position given by the virtual offset.
        # Note: pysam expects the virtual offset as is.
        inbam.seek(last_vo)
        with pysam.AlignmentFile(output_bam, "wb", header=inbam.header) as outbam:
            count = 0
            for read in inbam:
                # Filter to keep only unmapped reads.
                if read.is_unmapped:
                    outbam.write(read)
                    count += 1
            print(f"Extracted {count} unmapped reads to {output_bam}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract unmapped reads from a BAM file by seeking to the last mapped chunk using the BAI index."
    )
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("bai_file", help="Path to the corresponding BAI file")
    parser.add_argument("output_bam", help="Path for the output BAM file (unmapped reads)")
    args = parser.parse_args()

    extract_unmapped_reads_from_offset(args.bam_file, args.bai_file, args.output_bam)
