# extract-unmapped

A high-performance Python script to extract unmapped reads from a BAM file. This tool achieves its speed by making a critical assumption about the BAM file's structure, which allows it to skip reading the majority of the file.

**Note: This tool is optimized for speed by making assumptions about BAM file structure. Please review the usage guidelines to ensure it fits your use case.**

## Core Idea

The central concept is to avoid a full, linear scan of the entire BAM file. The script is built on the assumption that all unmapped reads are physically located at the end of the BAM file, after the last block containing a mapped read.

By making this assumption, the script can use the BAM Index (.bai) file to find the location of the last mapped read and start reading only from that point onwards, saving a significant amount of I/O and processing time on large files.

## How It Works (Implementation)

The script performs the following steps:

1. **Parse BAI Index**: It reads the binary .bai file to inspect the index data. The BAI contains pointers (virtual offsets) to compressed data blocks (BGZF blocks) in the BAM file that contain mapped reads.

2. **Find Max Offset**: It iterates through all references and all index bins to find the maximum `chunk_end` virtual offset. This offset points to the end of the BGZF block containing the last known mapped read in the entire file.

3. **Seek in BAM**: The script opens the main BAM file and performs a seek operation directly to this maximum virtual offset. This effectively jumps over all the mapped read data.

4. **Extract and Filter**: From that position, it reads until the end of the file. It inspects each alignment record and writes only the reads flagged as unmapped to a new output BAM file, preserving the original header.

## Limitations and Use Cases

This tool provides significant performance benefits for specific scenarios. Understanding when to use it ensures optimal results:

### Important Consideration: Paired-End Data

In coordinate-sorted BAM files from paired-end sequencing:

- Unmapped reads with mapped mates are typically stored adjacent to their mapped partner
- This means some unmapped reads may be interspersed throughout the file rather than concentrated at the end
- For such files, this tool will extract unmapped reads that appear after the last mapped read, but may miss those interspersed earlier

### File Structure Dependencies

The tool works best with files that maintain the original aligner output order:

- Direct output from aligners like BWA-MEM or Bowtie2 typically places unmapped reads at the end
- Files that have been merged, re-sorted, or post-processed may have different read distributions

## Usage Guidelines

### ✅ Optimal Use Cases

This tool provides excellent performance benefits for:

- Large BAM files with unmapped reads concentrated at the end (common with many aligners)
- Single-end sequencing data where unmapped reads are typically grouped together
- Files directly from aligners like BWA-MEM or Bowtie2 in their original output order
- Scenarios where speed is prioritized and you can verify completeness if needed

### ⚠️ Consider Alternatives For

- Paired-end data where complete extraction of all unmapped reads is critical
- Files that have been merged or extensively post-processed
- Workflows requiring guaranteed capture of every unmapped read regardless of file structure

## Installation and Dependencies

The script requires the `pysam` library.

```bash
pip install pysam
```

## Usage

```bash
python extract_unmapped_from_offset.py <input.bam> <input.bam.bai> <output_unmapped.bam>
```

- `input.bam`: Path to the input BAM file.
- `input.bam.bai`: Path to the corresponding BAM index file.
- `output_unmapped.bam`: Path for the new BAM file that will contain the extracted unmapped reads.

## Alternative: Complete Extraction

For guaranteed extraction of all unmapped reads regardless of file structure, use the standard samtools approach:

```bash
# Extract all unmapped reads with a full file scan
samtools view -f 4 -b -o unmapped.bam original.bam
```

This method scans the entire file and captures every unmapped read, trading speed for completeness.