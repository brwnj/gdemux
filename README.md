# Grouped Demultiplexing

We need to demultiplex, but want groups of barcodes to be joined into a single FASTQ file. And we want it to be easy.

This will demultiplex FASTQs using `fastq-multx` (`conda install -c bioconda fastq-multx`) then `cat` them into grouped
FASTQs rather than individual samples. The grouped FASTQs are
validated by total count against the individual sample sum.
This effectively subsets 'Undetermined' into smaller
groups of 'Undetermined' files.

# Group Named Output Files

Setting `--output-action` to "groupid" and running:

```
$ gdemux -a groupid -o out test_R1.fastq test_barcodes.txt
[2016-08-03 17:00 INFO] Found 10 samples across 3 groups within test_barcodes.txt
[2016-08-03 17:00 INFO] Demultiplexing (mismatches=1, distance=2, quality=0)
[2016-08-03 17:00 INFO] Joining reads across groups
[2016-08-03 17:00 INFO] Validating group read counts with sample counts
[2016-08-03 17:00 INFO] Processing complete

$ tree out
out
├── group1_I1.fastq
├── group1_R1.fastq
├── group1_R2.fastq
├── group2_I1.fastq
├── group2_R1.fastq
├── group2_R2.fastq
├── group3_I1.fastq
├── group3_R1.fastq
└── group3_R2.fastq
```

# Undetermined Subsets as Output Files

Setting `--output-action` to "undetermined" and running:

```
$ gdemux -a undetermined -o out test_R1.fastq test_barcodes.txt
[2016-08-03 17:01 INFO] Found 10 samples across 3 groups within test_barcodes.txt
[2016-08-03 17:01 INFO] Demultiplexing (mismatches=1, distance=2, quality=0)
[2016-08-03 17:01 INFO] Joining reads across groups
[2016-08-03 17:01 INFO] Validating group read counts with sample counts
[2016-08-03 17:01 INFO] Processing complete

$ tree out
out
├── group1
│   ├── Undetermined_I1.fastq
│   ├── Undetermined_R1.fastq
│   └── Undetermined_R2.fastq
├── group2
│   ├── Undetermined_I1.fastq
│   ├── Undetermined_R1.fastq
│   └── Undetermined_R2.fastq
└── group3
    ├── Undetermined_I1.fastq
    ├── Undetermined_R1.fastq
    └── Undetermined_R2.fastq
```

# Example Metadata

| groupid | sampleid | barcode      |
|---------|----------|--------------|
| group1  | McP-F1   | AAGGCGCTCCTT |
| group1  | McP-F2   | GATCTAATCGAG |
| group1  | McP-F3   | CTGATGTACACG |
| group2  | McP-F4   | ACGTATTCGAAG |
| group2  | CoP-F1   | GACGTTAAGAAT |
| group2  | CoP-F2   | TGGTGGAGTTTC |
| group3  | CoP-F3   | TTAACAAGGCAA |
| group3  | CoP-F4   | AACCGCATAAGT |
| group3  | McA-F1   | CCACAACGATCA |
| group3  | McA-F2   | AGTTCTCATTAA |

Extra columns can exist and differing column names can be used though they will need to be specified on the command line as `--sample-id`, `--group-id`, and `--barcode`.

A header isn't necessary either, though you'll need to specify more options. `--no-header` will be necessary, along with 0-based integers for the 3 columns, e.g. `--no-header --group-id 0 --sample-id 1 --barcode 2`.
