# Synthetic example data

Written by `make_examples.py`; every value is invented.

- `demo.vcf`: five made-up samples (sample IDs contain an underscore, which is what `--double-id` protects) and two variants.
- `demo.bed/.bim/.fam`: the same genotypes as a PLINK binary fileset.
- `truncated.*`: the same fileset with the last `.bed` byte removed, so validation fails.
