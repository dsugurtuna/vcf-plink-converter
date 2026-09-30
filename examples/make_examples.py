"""Write the synthetic example files in this folder.

python examples/make_examples.py

Every value is invented. demo.bed follows the PLINK 1 SNP-major format: three
magic bytes, then for each variant the samples packed four to a byte, two bits
each, first sample in the lowest bits (00 hom A1, 01 missing, 10 het, 11 hom
A2).
"""

from pathlib import Path

HERE = Path(__file__).parent
SAMPLES = ["SYN_01", "SYN_02", "SYN_03", "SYN_04", "SYN_05"]
VARIANTS = [("1", "rs900001", "10500", "G", "A"), ("1", "rs900002", "20750", "T", "C")]
# Genotypes as ALT-allele counts (None = missing), one row per variant.
GENOTYPES = [[0, 1, 2, 1, None], [1, 0, 0, 2, 1]]
VCF_GT = {0: "0/0", 1: "0/1", 2: "1/1", None: "./."}
# In the .bim, A1 = ALT and A2 = REF (as PLINK writes with --keep-allele-order).
BED_CODE = {2: 0b00, None: 0b01, 1: 0b10, 0: 0b11}


def bed_bytes(genotypes: list[list[int | None]]) -> bytes:
    out = bytearray(b"\x6c\x1b\x01")
    for row in genotypes:
        for start in range(0, len(row), 4):
            byte = 0
            for i, gt in enumerate(row[start : start + 4]):
                byte |= BED_CODE[gt] << (2 * i)
            out.append(byte)
    return bytes(out)


def main() -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=248956422>",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + SAMPLES
        ),
    ]
    for (chrom, vid, pos, ref, alt), row in zip(VARIANTS, GENOTYPES, strict=True):
        calls = [VCF_GT[g] for g in row]
        lines.append(
            "\t".join([chrom, vid, pos, ref, alt, ".", "PASS", ".", "GT", *calls])
        )
    (HERE / "demo.vcf").write_text("\n".join(lines) + "\n")

    fam = "".join(f"{s} {s} 0 0 0 -9\n" for s in SAMPLES)
    bim = "".join(f"{c}\t{v}\t0\t{p}\t{alt}\t{ref}\n" for c, v, p, ref, alt in VARIANTS)
    for name, bed in (
        ("demo", bed_bytes(GENOTYPES)),
        ("truncated", bed_bytes(GENOTYPES)[:-1]),
    ):
        (HERE / f"{name}.fam").write_text(fam)
        (HERE / f"{name}.bim").write_text(bim)
        (HERE / f"{name}.bed").write_bytes(bed)


if __name__ == "__main__":
    main()
