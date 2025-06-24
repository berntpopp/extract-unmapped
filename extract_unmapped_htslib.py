# extract_unmapped_htslib.py
# Extract unmapped reads from a BAM file using direct htslib C library calls via ctypes.
# This script replicates the functionality of extract_unmapped_from_offset.py without pysam dependency.

import argparse
import ctypes
import ctypes.util
import os
import sys


# Load htslib shared library
def load_htslib():
    """Load the htslib shared library using ctypes."""
    lib_name = ctypes.util.find_library("hts")
    if lib_name is None:
        raise RuntimeError(
            "htslib shared library not found. Please install htslib:\n"
            "  Ubuntu/Debian: apt-get install libhts-dev\n"
            "  CentOS/RHEL: yum install htslib-devel\n"
            "  macOS: brew install htslib"
        )
    try:
        return ctypes.CDLL(lib_name)
    except OSError as e:
        raise RuntimeError(f"Failed to load htslib library '{lib_name}': {e}")

# Load the library
hts = load_htslib()


# Define C data structures
class SamFile(ctypes.Structure):
    """samFile structure with BGZF pointer access"""
    _fields_ = [
        ("fp", ctypes.c_void_p),  # BGZF pointer for seeking
        # Other fields exist but we only need fp for our seeking operation
    ]

class BamHdr(ctypes.Structure):
    """Opaque structure for BAM header"""
    pass

class BGZF(ctypes.Structure):
    """Opaque structure for BGZF handle"""
    pass

class Bam1Core(ctypes.Structure):
    """BAM alignment core structure - we only need fields up to flag"""
    _fields_ = [
        ("tid", ctypes.c_int32),        # chromosome ID
        ("pos", ctypes.c_int32),        # 0-based leftmost coordinate
        ("bin", ctypes.c_uint16),       # bin calculated by bam_reg2bin()
        ("qual", ctypes.c_uint8),       # mapping quality
        ("l_qname", ctypes.c_uint8),    # length of the query name
        ("flag", ctypes.c_uint16),      # bitwise flag
        ("n_cigar", ctypes.c_uint16),   # number of CIGAR operations
        ("l_qseq", ctypes.c_int32),     # length of the query sequence
        ("mtid", ctypes.c_int32),       # chromosome ID of next read in template
        ("mpos", ctypes.c_int32),       # 0-based leftmost coordinate of next read
        ("isize", ctypes.c_int32),      # observed template length
    ]

class Bam1(ctypes.Structure):
    """BAM alignment structure"""
    _fields_ = [
        ("core", Bam1Core),             # core information
        ("l_data", ctypes.c_int),       # current length of bam1_t::data
        ("m_data", ctypes.c_int),       # maximum length of bam1_t::data
        ("data", ctypes.POINTER(ctypes.c_uint8)),  # all variable-length data
    ]

# Create pointer types
SamFile_p = ctypes.POINTER(SamFile)
BamHdr_p = ctypes.POINTER(BamHdr)
BGZF_p = ctypes.POINTER(BGZF)
Bam1_p = ctypes.POINTER(Bam1)

# Define C function prototypes
# File operations
hts.hts_open.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
hts.hts_open.restype = SamFile_p

hts.hts_close.argtypes = [SamFile_p]
hts.hts_close.restype = ctypes.c_int

# Header operations
hts.sam_hdr_read.argtypes = [SamFile_p]
hts.sam_hdr_read.restype = BamHdr_p

hts.sam_hdr_write.argtypes = [SamFile_p, BamHdr_p]
hts.sam_hdr_write.restype = ctypes.c_int

hts.sam_hdr_destroy.argtypes = [BamHdr_p]
hts.sam_hdr_destroy.restype = None

# BAM record operations
hts.bam_init1.argtypes = []
hts.bam_init1.restype = Bam1_p

hts.bam_destroy1.argtypes = [Bam1_p]
hts.bam_destroy1.restype = None

hts.sam_read1.argtypes = [SamFile_p, BamHdr_p, Bam1_p]
hts.sam_read1.restype = ctypes.c_int

hts.sam_write1.argtypes = [SamFile_p, BamHdr_p, Bam1_p]
hts.sam_write1.restype = ctypes.c_int

# SamFile structure already defined above

# BGZF seeking functions
hts.bgzf_seek.argtypes = [ctypes.c_void_p, ctypes.c_int64, ctypes.c_int]
hts.bgzf_seek.restype = ctypes.c_int

# Define BAM constants
BAM_FUNMAP = 4  # read unmapped flag


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


def extract_unmapped_with_htslib(bam_file, bai_file, output_bam, debug=False):
    """
    Extract unmapped reads from a BAM file using direct htslib calls.
    Uses BAI index to find the maximum virtual offset and seeks to that position.
    """
    # Get the maximum virtual offset from the BAI file
    last_vo = get_last_chunk_end(bai_file)
    print(f"Last mapped virtual offset (from BAI): {last_vo}")
    
    if debug:
        print(f"DEBUG: BAI file size: {os.path.getsize(bai_file)} bytes")
        print(f"DEBUG: BAM file size: {os.path.getsize(bam_file)} bytes")

    # Initialize all resources to None for cleanup
    in_fp = None
    out_fp = None
    header = None
    record = None

    try:
        if debug:
            print("DEBUG: Opening input BAM file...")
        
        # Open the input BAM file
        in_fp = hts.hts_open(bam_file.encode(), b"r")
        if not in_fp:
            raise IOError(f"Failed to open input BAM file: {bam_file}")
        
        if debug:
            print("DEBUG: Input BAM file opened successfully")

        # Read the header
        if debug:
            print("DEBUG: Reading BAM header...")
        header = hts.sam_hdr_read(in_fp)
        if not header:
            raise IOError(f"Failed to read header from BAM file: {bam_file}")
        
        if debug:
            print("DEBUG: BAM header read successfully")

        # Seek to the position given by the virtual offset
        if debug:
            print(f"DEBUG: Attempting to seek to virtual offset {last_vo}")
        
        # Access the BGZF pointer from the samFile structure
        try:
            bgzf_ptr = in_fp.contents.fp
            if debug:
                print(f"DEBUG: Got BGZF pointer: {bgzf_ptr}")
            
            seek_result = hts.bgzf_seek(bgzf_ptr, last_vo, 0)  # SEEK_SET = 0
            if seek_result < 0:
                raise IOError(f"Failed to seek to virtual offset {last_vo}")
            
            if debug:
                print(f"DEBUG: Seek successful, result: {seek_result}")
        except Exception as e:
            if debug:
                print(f"DEBUG: Seek failed with error: {e}")
            raise IOError(f"Failed to seek to virtual offset {last_vo}: {e}")

        # Open the output BAM file
        if debug:
            print("DEBUG: Opening output BAM file...")
        out_fp = hts.hts_open(output_bam.encode(), b"wb")
        if not out_fp:
            raise IOError(f"Failed to open output BAM file: {output_bam}")
        
        if debug:
            print("DEBUG: Output BAM file opened successfully")

        # Write the header to the output file
        if debug:
            print("DEBUG: Writing header to output file...")
        if hts.sam_hdr_write(out_fp, header) < 0:
            raise IOError("Failed to write header to output BAM file")
        
        if debug:
            print("DEBUG: Header written successfully")

        # Initialize a reusable record object
        if debug:
            print("DEBUG: Initializing BAM record...")
        record = hts.bam_init1()
        if not record:
            raise IOError("Failed to initialize BAM record")
        
        if debug:
            print("DEBUG: BAM record initialized successfully")

        # Read and filter records
        count = 0
        read_count = 0
        if debug:
            print("DEBUG: Starting to read records...")
        
        while True:
            read_result = hts.sam_read1(in_fp, header, record)
            if read_result < 0:
                if debug:
                    print(f"DEBUG: End of file reached, sam_read1 returned {read_result}")
                break
            
            read_count += 1
            if debug and read_count <= 5:
                try:
                    flag = record.contents.core.flag
                    is_unmapped = (flag & BAM_FUNMAP) != 0
                    print(f"DEBUG: Read {read_count}: flag={flag}, unmapped={is_unmapped}")
                except Exception as e:
                    print(f"DEBUG: Error accessing record {read_count}: {e}")
                    break
            
            # Check if the read is unmapped
            try:
                if record.contents.core.flag & BAM_FUNMAP:
                    write_result = hts.sam_write1(out_fp, header, record)
                    if write_result < 0:
                        raise IOError("Failed to write record to output BAM file")
                    count += 1
                    if debug and count <= 3:
                        print(f"DEBUG: Wrote unmapped read {count}")
            except Exception as e:
                if debug:
                    print(f"DEBUG: Error processing record {read_count}: {e}")
                break

        if debug:
            print(f"DEBUG: Processed {read_count} total reads")
        print(f"Extracted {count} unmapped reads to {output_bam}")

    finally:
        # Clean up all resources
        if record:
            hts.bam_destroy1(record)
        if header:
            hts.sam_hdr_destroy(header)
        if in_fp:
            hts.hts_close(in_fp)
        if out_fp:
            hts.hts_close(out_fp)


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Extract unmapped reads from BAM files using direct htslib calls",
        epilog="""
Examples:
  %(prog)s input.bam input.bam.bai unmapped.bam

Performance Notes:
  This is an alternative implementation that uses direct htslib C library calls
  instead of pysam. It provides the same optimization as the pysam version but
  without the pysam dependency.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("bam_file", 
                       help="Path to the input BAM file")
    parser.add_argument("bai_file", 
                       help="Path to the corresponding BAI index file")
    parser.add_argument("output_bam", 
                       help="Path for the output BAM file (unmapped reads)")
    parser.add_argument("-d", "--debug", 
                       action="store_true",
                       help="Enable debug output")
    parser.add_argument("-v", "--version", 
                       action="version", 
                       version="%(prog)s 1.0 (htslib)")
    
    args = parser.parse_args()
    
    try:
        # Validate input files exist
        if not os.path.isfile(args.bam_file):
            raise IOError(f"Input BAM file not found: {args.bam_file}")
        if not os.path.isfile(args.bai_file):
            raise IOError(f"BAI index file not found: {args.bai_file}")
            
        extract_unmapped_with_htslib(args.bam_file, args.bai_file, args.output_bam, debug=args.debug)
        
    except (IOError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()