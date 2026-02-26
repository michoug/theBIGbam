# Preprocessing your data

## Mapping

You can map reads to a reference assembly using `mapping-per-sample`. The command generates a sorted .bam file along with its corresponding index (.bam.bai).

Available mapper choices: `minimap2-sr`, `bwa-mem2`, `minimap2-ont`, `minimap2-pb`, `minimap2-hifi`, `minimap2-no-preset`, `minimap2-sr-secondary`

The mapping strategy used here supports circularized mapping via the `--circular` flag. In this mode, each contig is concatenated to itself during mapping so that reads spanning the contig boundaries map seamlessly without being split at the ends. After mapping, the BAM is automatically converted to SAM-spec circular coordinates: `@SQ` headers carry the original (real) contig lengths, and read positions are normalised so that POS is always within the contig bounds. Ghost secondary and supplementary alignments created by the doubled reference are filtered out (the number of discarded ghost alignments is printed to stdout). If you use this option, ensure the `--circular` flag is also specified during database computation to correctly interpret the mappings.

### SAM specification for circular genomes

The output BAM follows the SAM specification (section 1.4) for circular references:

> *"POS plus the sum of the lengths of M/=/X/D/N CIGAR operations may exceed LN. Coordinates greater than LN are interpreted by subtracting LN so that bases at LN+1, LN+2, LN+3, … are considered to be mapped at positions 1, 2, 3, …; thus each (1-based) position p is interpreted as ((p − 1) mod LN) + 1."*

This means reads that cross the origin of a circular genome will have their alignment extend past the contig length in the BAM file, which is valid per the specification.

**Compatibility warning:** although this encoding is SAM-spec-compliant, most downstream tools (IGV, samtools depth, variant callers, etc.) do not yet implement circular coordinate handling. Programs that are unaware of this part of the specification will likely produce incorrect results for reads wrapping around the origin. Verify that your downstream tool explicitly supports circular genomes before relying on the output outside of theBIGbam.

**MAPQ caveat:** in circular mode the mapper sees the doubled reference as two distinct loci, so every read produces at least one ghost alignment that competes with the real one. The mapper therefore assigns lower MAPQ scores than the same read would receive against a single-copy reference, because MAPQ quantifies the confidence that the *reported* alignment is the only good one — and in a doubled reference there are always at least two.

The circular-BAM conversion step removes these ghost secondary and supplementary alignments, restoring the correct read count. For reads whose *only* competing alignments were ghosts (i.e. after removal they have no remaining secondary or supplementary records), the conversion sets MAPQ to 60, which is the conventional value for a uniquely mapping read. Reads that retain legitimate secondary or supplementary alignments after ghost removal keep their original MAPQ unchanged.

Full MAPQ recalculation is not possible because both minimap2 and bwa-mem2 derive MAPQ from internal values that are never written to BAM output:

- **minimap2** uses chaining scores (*f*₁, *f*₂) and anchor counts to compute MAPQ; these are discarded after alignment.
- **bwa-mem2** uses chaining scores, seed coverage, and the repetitive fraction (*frac_rep*); these are likewise internal only.

Without access to those intermediate values, no post-hoc tool can faithfully reproduce the mapper's MAPQ formula. Setting MAPQ=60 for de-ghosted unique mappers is the most accurate correction available.

Another feature of this mapping process is the computation of MD tags using samtools. MD tags are useful to quickly identify mismatches among the mapped reads.

You can pass extra parameters to the underlying mapper with `--minimap2-params` or `--bwa-params`. Extra parameters are appended after the preset flags, so if you specify a flag that already appears in the preset, your value overrides it (last-flag-wins behaviour in both minimap2 and bwa-mem2).

Examples:

```sh
DIR="examples/inputs/HK97"
thebigbam mapping-per-sample --mapper minimap2-sr \
  -r1 "${DIR}/HK97_R1_illumina.fastq.gz" \
  -r2 "${DIR}/HK97_R2_illumina.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -o "${DIR}/HK97_GCF_000848825.1_with_MD.bam" --circular

# Using interleaved paired-end FASTQ
thebigbam mapping-per-sample --mapper bwa-mem2 \
  --interleaved "${DIR}/HK97_interleaved.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -o "${DIR}/HK97_bwa.bam"

# Passing extra mapper parameters
thebigbam mapping-per-sample --mapper minimap2-sr \
  -r1 "${DIR}/HK97_R1_illumina.fastq.gz" \
  -r2 "${DIR}/HK97_R2_illumina.fastq.gz" \
  -a "${DIR}/HK97_GCF_000848825.1.fasta" \
  -o "${DIR}/HK97_custom.bam" \
  --minimap2-params "--secondary=no -N5"
```
