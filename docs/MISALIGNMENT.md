# Codon Change Annotation for Mismatches

When an annotation file is provided, each mismatch position inside a CDS is annotated with its codon impact: whether the base substitution is **synonymous** (same amino acid) or **non-synonymous** (different amino acid). Positions outside any CDS are labelled **Intergenic**.

## Per-read codon awareness

A key design choice is that codon changes are evaluated **per read**, not per position independently. This matters when a single read carries multiple mismatches within the same codon.

**Example:** Reference codon is `ATG` (Met). A read has mismatches at both position 1 (A→C) and position 2 (T→G). Evaluating each position alone would give mutant codons `CTG` and `AGG`, but the actual codon in that read is `CGG` (Arg). The per-read approach captures this correctly.

### Algorithm

For each read, during the MD tag walk:

1. **Collect mismatches**: all `(genome_position, read_base)` pairs from the MD tag are recorded.
2. **Map to CDS**: each mismatch position is looked up against a sorted index of CDS intervals (binary search). Positions not in any CDS are recorded as Intergenic.
3. **Group by codon**: mismatches falling in the same CDS and same codon (determined by `offset_in_cds / 3`) are grouped together.
4. **Build the mutant codon**: starting from the reference codon (extracted from the CDS nucleotide sequence), **all** mismatches from this read are substituted simultaneously. For reverse-strand CDS, bases are complemented before substitution.
5. **Translate and classify**: both reference and mutant codons are translated using the standard genetic code. If the amino acid is unchanged, the change is Synonymous; otherwise Non-synonymous.
6. **Record**: the resulting `(category, mutant_codon, amino_acid_change)` tuple is counted at each affected genome position.

### Dominant selection

After processing all reads for a contig, each genome position may have accumulated multiple distinct codon change variants (from different read subpopulations). The variant with the **highest read count** is selected as the dominant one and stored in the database.

**Example with subpopulations:** Consider a codon where position B is mutated in all variant reads, but position A is only mutated in some of them:

- 50 reads: only B mutated → mutant codon `TYG`, translates to Val (Synonymous)
- 30 reads: A and B mutated → mutant codon `XYG`, translates to Phe (Non-synonymous)

At position B, both variants compete. The B-alone variant wins (50 > 30), so position B is stored as Synonymous with codon `TYG`. Position A only has the double-mutant variant (30 reads), so it is stored as Non-synonymous with codon `XYG`.

### Two independent dominant selections

Each mismatch position stores two kinds of dominant information that are computed independently:

- **`Sequence` / `Sequence_prevalence`** (tooltip: Sequence / Prevalence): the dominant **nucleotide** at this position, selected by counting how many reads carry each alternative base. This is a per-position metric — it only looks at the single nucleotide change at this position, regardless of what happens at neighbouring positions in the same codon.

- **`Codon_change` / `AA_change`** (tooltip: Codon / Amino acid): the dominant **codon variant**, selected by counting how many reads carry each distinct whole-codon combination. This is a per-codon metric — it considers all mismatches within the codon simultaneously.

In nearly all cases the two agree: the dominant nucleotide at a position will be the same base found in the dominant codon. They can differ in rare situations where read subpopulations carry different combinations of mutations within the same codon. In the example above, all 80 reads agree on the same base change at position B (so `Sequence` reflects 80 reads), but only 50 of those reads share the same codon context (so `Codon_change` reflects the 50-read subpopulation). The nucleotide in the dominant codon will still match `Sequence` here, but in more complex scenarios with three or more subpopulations they could theoretically disagree.

## Database columns

The following columns are added to the `Feature_mismatches` table:

| Column                | Tooltip label | Description                                                                |
| --------------------- | ------------- | -------------------------------------------------------------------------- |
| `Sequence`            | Sequence      | Dominant alternative nucleotide at this position                           |
| `Sequence_prevalence` | Prevalence    | Percentage of reads carrying that nucleotide (relative to coverage)        |
| `Codon_category`      | Category      | `Synonymous`, `Non-synonymous`, or `Intergenic`                            |
| `Codon_change`        | Codon         | The dominant mutant codon (e.g. `ACG`). NULL for Intergenic positions      |
| `AA_change`           | Amino acid    | The resulting amino acid, e.g. `V (Valine)`. NULL for Intergenic positions |

Codon columns are NULL for mismatch positions where no dominant base was identified (i.e. positions stored without a `Sequence` value).

## CDS index

CDS intervals are extracted from the annotation file during parsing. For each contig, a sorted list of CDS intervals (start, end, strand, nucleotide sequence) is built. Lookup uses binary search on CDS start positions, scanning backwards to handle overlapping/nested CDS features correctly. Only longest isoforms are included when locus tag information is available.
