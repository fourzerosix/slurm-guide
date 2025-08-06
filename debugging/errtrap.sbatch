#!/bin/bash
#
# Author: Dolphin Whisperer
# Email: jeremy.bell@nih.gov
# Created: 2025-07-18
# Description: This script exemplifies catching any unexpected error and report exactly where it occurred,
# combining set -eEuo pipefail with an ERR trap - it builds an index, runs an alignment, produces a pretty output file and makes it home before curfew. 
# set -Eeuo pipefail:
#     -e: exit on error
#     -E: propagate ERR traps into subshells/functions
#     -u: error on unset variables
#      pipefail: fail on pipeline error
# The trap prints the script name, line number (top of the faulty frame), and the failed command—then exits.
# Clear, automatic failure reporting with context—no guessing which command failed.

#SBATCH --job-name=issa-trap
#SBATCH --output=issa-trap_%j.out
#SBATCH --error=issa-trap_%j.err

set -Eeuo pipefail
trap 'echo "[$(date)] ERROR in ${BASH_SOURCE[0]} at line ${BASH_LINENO[0]}: command \"${BASH_COMMAND}\" failed."; exit 1' ERR

# Load necessary modules
ml bwa
ml samtools

# Create dummy reference genome
cat <<EOF > ref.fa
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF

# Create dummy FASTQ reads
cat <<EOF > reads.fq
@read1
ACGTACGTACGT
+
FFFFFFFFFFFF
@read2
ACGTACGTACGT
+
FFFFFFFFFFFF
EOF

# Build BWA index (required for bwa mem)
bwa index ref.fa

# Run alignment and sort the result
bwa mem ref.fa reads.fq > aln.sam
samtools sort -O BAM aln.sam -o aln.sorted.bam
