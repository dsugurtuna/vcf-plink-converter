"""Malformed inputs must be rejected by FileValidator and VCFInspector."""

import gzip
from pathlib import Path

from vcf_converter.inspector import VCFInspector
from vcf_converter.validator import FileValidator

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS_1\tS_2\n"


def _fileset(tmp_path: Path, bed: bytes, n_variants: int = 1) -> str:
    prefix = tmp_path / "x"
    Path(f"{prefix}.bed").write_bytes(bed)
    Path(f"{prefix}.bim").write_text(
        "".join(f"1\trs{i}\t0\t{i}00\tA\tG\n" for i in range(n_variants))
    )
    Path(f"{prefix}.fam").write_text("F1 I1 0 0 1 -9\nF2 I2 0 0 2 -9\n")
    return str(prefix)


def test_truncated_bed_is_rejected(tmp_path: Path) -> None:
    # 3 variants x 2 samples needs 3 + 3 * 1 = 6 bytes.
    prefix = _fileset(tmp_path, b"\x6c\x1b\x01\x00\x00", n_variants=3)
    report = FileValidator().validate_plink_binary(prefix)
    assert not report.all_valid
    assert "needs 6" in report.invalid[0]


def test_consistent_fileset_passes(tmp_path: Path) -> None:
    prefix = _fileset(tmp_path, b"\x6c\x1b\x01\x00\x00\x00", n_variants=3)
    assert FileValidator().validate_plink_binary(prefix).valid == 3


def test_bim_with_wrong_columns_is_rejected(tmp_path: Path) -> None:
    prefix = _fileset(tmp_path, b"\x6c\x1b\x01\x00")
    Path(f"{prefix}.bim").write_text("1 rs1 100 A G\n")
    report = FileValidator().validate_plink_binary(prefix)
    assert any("without 6 columns" in m for m in report.invalid)


def test_vcf_without_mandatory_columns_is_rejected(tmp_path: Path) -> None:
    vcf = tmp_path / "x.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\n1\t100\trs1\n")
    report = FileValidator().validate_vcf(vcf)
    assert "8 mandatory columns" in report.invalid[0]


def test_gzipped_vcf_and_format_id_without_comma(tmp_path: Path) -> None:
    vcf = tmp_path / "x.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT>\n##contig=<ID=chr1>\n")
        fh.write(HEADER + "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n")
    assert FileValidator().validate_vcf(vcf).all_valid
    info = VCFInspector().inspect(vcf)
    assert info.format_fields == ["GT"]
    assert info.contigs == ["chr1"]
    assert info.sample_names == ["S_1", "S_2"]
