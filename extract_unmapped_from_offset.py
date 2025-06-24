#!/usr/bin/env python3
"""
extract_unmapped_from_offset.py

High-performance extraction of unmapped reads from BAM files using BAI index optimization.

This script uses the BAM index (BAI) to find the maximum virtual offset among all mapped 
regions, then seeks directly to that position to extract unmapped reads that appear after 
the last mapped read. This provides significant performance benefits for large files where 
unmapped reads are concentrated at the end.

Performance Trade-off:
- Speed: Much faster than full file scans for appropriate use cases
- Completeness: May miss unmapped reads interspersed with mapped reads in paired-end data

Dependencies: pysam
"""

import argparse
import os
import sys
import pysam

def read_uint32(f):
    """Read a 32-bit unsigned integer from binary file in little-endian format."""
    return int.from_bytes(f.read(4), byteorder="little", signed=False)

def read_uint64(f):
    """Read a 64-bit unsigned integer from binary file in little-endian format."""
    return int.from_bytes(f.read(8), byteorder="little", signed=False)

def get_last_chunk_end(bai_filename):
    """
    Parse BAI index file to find the maximum virtual offset among all mapped regions.
    
    Args:
        bai_filename (str): Path to the BAM index (.bai) file
        
    Returns:
        int: Maximum virtual offset (chunk_end) found in the index
        
    Raises:
        IOError: If the BAI file cannot be read or is malformed
    """
    max_vo = 0
    try:
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
    except (IOError, OSError) as e:
        raise IOError(f"Failed to read BAI file '{bai_filename}': {e}")
    return max_vo

def extract_unmapped_reads_from_offset(bam_file, bai_file, output_bam, verbose=True, debug=False):
    """
    Extract unmapped reads from a BAM file using BAI index optimization.
    
    This function uses the BAI index to find the last mapped read position,
    then seeks to that position and extracts all subsequent unmapped reads.
    
    Args:
        bam_file (str): Path to the input BAM file
        bai_file (str): Path to the corresponding BAI index file
        output_bam (str): Path for the output BAM file containing unmapped reads
        verbose (bool): Whether to print progress information
        
    Returns:
        int: Number of unmapped reads extracted
        
    Raises:
        IOError: If any file operations fail
        ImportError: If pysam is not available
    """
    # Validate input files exist
    if not os.path.isfile(bam_file):
        raise IOError(f"Input BAM file not found: {bam_file}")
    if not os.path.isfile(bai_file):
        raise IOError(f"BAI index file not found: {bai_file}")
    
    # Get the maximum virtual offset from the BAI file
    last_vo = get_last_chunk_end(bai_file)
    if verbose:
        print(f"Last mapped virtual offset (from BAI): {last_vo}")
    if debug:
        print(f"DEBUG: BAI file size: {os.path.getsize(bai_file)} bytes")
        print(f"DEBUG: BAM file size: {os.path.getsize(bam_file)} bytes")

    try:
        # Open the BAM file and seek to the computed virtual offset
        with pysam.AlignmentFile(bam_file, "rb") as inbam:
            if debug:
                print(f"DEBUG: Opened BAM file successfully")
                print(f"DEBUG: BAM header references: {inbam.nreferences}")
            
            # Seek to the position given by the virtual offset
            # Note: pysam expects the virtual offset as is
            inbam.seek(last_vo)
            if debug:
                print(f"DEBUG: Seeked to virtual offset {last_vo}")
            
            with pysam.AlignmentFile(output_bam, "wb", header=inbam.header) as outbam:
                if debug:
                    print(f"DEBUG: Created output BAM file")
                
                count = 0
                read_count = 0
                for read in inbam:
                    read_count += 1
                    if debug and read_count <= 5:
                        print(f"DEBUG: Read {read_count}: mapped={not read.is_unmapped}, flag={read.flag}")
                    
                    # Filter to keep only unmapped reads
                    if read.is_unmapped:
                        outbam.write(read)
                        count += 1
                        if debug and count <= 3:
                            print(f"DEBUG: Wrote unmapped read {count}")
                
                if debug:
                    print(f"DEBUG: Processed {read_count} total reads")
                if verbose:
                    print(f"Extracted {count} unmapped reads to {output_bam}")
                return count
                
    except Exception as e:
        raise IOError(f"Failed to process BAM files: {e}")

def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Extract unmapped reads from BAM files using BAI index optimization",
        epilog="""
Examples:
  %(prog)s input.bam input.bam.bai unmapped.bam
  %(prog)s --quiet large_file.bam large_file.bam.bai output/unmapped.bam

Performance Notes:
  This tool is optimized for speed by assuming unmapped reads are located
  after the last mapped read. It may miss unmapped reads interspersed with
  mapped reads in paired-end data. For complete extraction, use:
  samtools view -f 4 -b input.bam > unmapped.bam
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("bam_file", 
                       help="Path to the input BAM file")
    parser.add_argument("bai_file", 
                       help="Path to the corresponding BAI index file")
    parser.add_argument("output_bam", 
                       help="Path for the output BAM file (unmapped reads)")
    parser.add_argument("-q", "--quiet", 
                       action="store_true",
                       help="Suppress progress output")
    parser.add_argument("-d", "--debug", 
                       action="store_true",
                       help="Enable debug output")
    parser.add_argument("-v", "--version", 
                       action="version", 
                       version="%(prog)s 1.0")
    
    args = parser.parse_args()
    
    try:
        count = extract_unmapped_reads_from_offset(
            args.bam_file, 
            args.bai_file, 
            args.output_bam, 
            verbose=not args.quiet,
            debug=args.debug
        )
        if not args.quiet:
            print(f"Success: Extracted {count} unmapped reads")
    except (IOError, ImportError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
