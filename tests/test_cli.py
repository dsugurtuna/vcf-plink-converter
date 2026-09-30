"""CLI tests against examples/ (no PLINK needed)."""

from pathlib import Path

import pytest

from vcf_converter.__main__ import main

EX = Path(__file__).resolve().parent.parent / "examples"


def test_inspect(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["inspect", str(EX / "demo.vcf")]) == 0
    out = capsys.readouterr().out
    assert "samples:  5 (SYN_01" in out
    assert "variants: 2" in out


def test_validate_good_and_truncated(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["validate", str(EX / "demo.vcf"), str(EX / "demo")]) == 0
    assert main(["validate", str(EX / "truncated")]) == 1
    assert ".bed is 6 bytes; 2 variants x 5 samples needs 7" in capsys.readouterr().out


def test_convert_without_plink_fails_cleanly(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(
        [
            "vcf-to-plink",
            str(EX / "demo.vcf"),
            str(tmp_path / "o"),
            "--plink",
            "no-plink",
        ]
    )
    assert code == 1
    assert capsys.readouterr().out.startswith("failed:")
