"""PLINK command lines and result handling, with PLINK replaced by a fake."""

import gzip
import subprocess
from pathlib import Path

from vcf_converter.converter import FormatConverter

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS_1\tS_2\n"


def test_vcf_import_keeps_whole_sample_ids() -> None:
    cmd = FormatConverter()._build_vcf_to_plink_cmd("in.vcf", "out")
    assert "--double-id" in cmd
    assert "--keep-allele-order" in cmd


def test_plink_to_vcf_gz_uses_bgz_and_counts(tmp_path: Path) -> None:
    def fake(argv: list[str]) -> subprocess.CompletedProcess[str]:
        out = argv[argv.index("--out") + 1]
        with gzip.open(f"{out}.vcf.gz", "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n" + HEADER)
            fh.write("1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n")
        return subprocess.CompletedProcess(argv, 0, "", "")

    result = FormatConverter(runner=fake).plink_to_vcf("in", str(tmp_path / "o.vcf.gz"))
    assert result.success
    assert result.output_path.endswith("o.vcf.gz")
    assert (result.sample_count, result.variant_count) == (2, 1)


def test_missing_plink_is_reported() -> None:
    result = FormatConverter(plink_binary="plink-does-not-exist").vcf_to_plink(
        "in.vcf", "out"
    )
    assert not result.success
    assert "plink-does-not-exist" in result.message
