# Why it's built this way

## The problem

Converting between VCF and PLINK binary files is routine, but PLINK 1.9's defaults change sample IDs and can change which allele is called REF, and a damaged `.bed` file is not noticed until something downstream fails. The conversion should give back what went in, and a bad file should be caught at the door.

## Design choices

**Why `--double-id` on import?** Because by default PLINK 1.9 splits a VCF sample ID at a single underscore into family and individual IDs, and rejects IDs with more than one. Running the example data through PLINK without the flag turns `SYN_01` into FID `SYN`, IID `01`. Many biobank sample IDs contain underscores.

**Why `vcf-iid` on export?** Because plain `--recode vcf` writes `FID_IID` as the sample name. With `--double-id` on import that doubles the ID (`SYN_01_SYN_01`); `vcf-iid` writes it once.

**Why `--keep-allele-order` both ways?** Because PLINK 1.9 may otherwise make the major allele A2, and on export A2 is written as REF. A REF that follows frequency rather than the reference genome breaks tools that check REF against the genome.

**Why check the `.bed` size?** Because it is fully determined: three magic bytes plus, for each variant, one byte per four samples. One comparison with the `.bim` and `.fam` row counts catches a truncated copy or a `.bed` paired with the wrong `.bim`, which the magic bytes alone never would.

**Why run PLINK through an injectable runner?** Because command lines and failure handling can then be tested anywhere. It is not a substitute for running the real tool: doing so found a bug (POS and ID swapped in the example VCF) that every structural check had passed.

**Why no dependencies for validation and inspection?** So they can run as a gate on any machine, including CI, without PLINK installed.

## Questions worth asking

**"Does `--keep-allele-order` guarantee REF matches the reference genome?"**
No. It preserves whatever REF the input had. If the input VCF came from a tool that already wrote correct REF alleles, it stays correct. If the data started as PLINK files, REF may be wrong, and the fix is PLINK's `--a2-allele` with a reference VCF (and `--real-ref-alleles`), which is on the roadmap.

**"What gets lost in the conversion?"**
Everything the PLINK 1 binary format cannot hold: phase, dosages or genotype probabilities, FORMAT fields other than GT, filters and INFO annotations, and multi-allelic sites. For imputed data you would keep the VCF (or use PLINK 2 formats) rather than round-trip through `.bed`.

**"Why not use PLINK 2?"**
PLINK 2 handles several of these issues better (for example, it keeps allele order by default). Much existing tooling and many pipelines still expect PLINK 1.9 behaviour, so this tool makes 1.9 behave safely rather than assuming everyone has moved.

## What's next

- An optional `--a2-allele` step to set REF from a reference VCF.
- Check chromosome naming (`1` versus `chr1`) against a target build.
- Record checksums of inputs and outputs for each conversion.
