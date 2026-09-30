"""Structural checks for VCF and PLINK binary files before conversion."""

from __future__ import annotations

import gzip
import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .inspector import iter_vcf_lines

VCF_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


@dataclass
class ValidationReport:
    """Validation outcome for a set of files."""

    files_checked: int = 0
    valid: int = 0
    invalid: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def all_valid(self) -> bool:
        return len(self.invalid) == 0


class FileValidator:
    """Validate VCF and PLINK binary files.

    VCF: the first line must declare ``##fileformat=VCF``, and a ``#CHROM``
    header with the eight mandatory columns must come before any data.
    PLINK: ``.bed``, ``.bim`` and ``.fam`` must exist, ``.bim``/``.fam`` rows
    need six columns, the ``.bed`` must start with the SNP-major magic bytes
    and its size must equal ``3 + variants * ceil(samples / 4)``, which
    catches truncated or mismatched filesets.
    """

    VCF_HEADER = "##fileformat=VCF"
    PLINK_BED_MAGIC = b"\x6c\x1b\x01"

    def validate_vcf(self, vcf_path: str | Path) -> ValidationReport:
        report = ValidationReport(files_checked=1)
        path = Path(vcf_path)
        if not path.exists():
            report.invalid.append(f"File not found: {vcf_path}")
            return report
        try:
            lines = iter_vcf_lines(path)
            first = next(lines, "")
            if not first.startswith(self.VCF_HEADER):
                report.invalid.append(f"Missing ##fileformat=VCF line: {vcf_path}")
                return report
            for line in lines:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    if line.split("\t")[:8] != VCF_COLUMNS:
                        report.invalid.append(
                            f"#CHROM header lacks the 8 mandatory columns: {vcf_path}"
                        )
                        return report
                    report.valid = 1
                    return report
                break
            report.invalid.append(f"No #CHROM header line before data: {vcf_path}")
        except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as exc:
            report.invalid.append(f"Read error: {exc}")
        return report

    @staticmethod
    def _rows_with_columns(path: Path, expected: int) -> tuple[int, int]:
        """Return (rows, rows with the wrong number of columns)."""
        rows = bad = 0
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                rows += 1
                if len(line.split()) != expected:
                    bad += 1
        return rows, bad

    def validate_plink_binary(self, bfile_prefix: str | Path) -> ValidationReport:
        bed, bim, fam = (
            Path(f"{bfile_prefix}{ext}") for ext in (".bed", ".bim", ".fam")
        )
        report = ValidationReport(files_checked=3)
        for fp, label in ((bed, "bed"), (bim, "bim"), (fam, "fam")):
            if not fp.exists():
                report.invalid.append(f"Missing .{label} file: {fp}")

        counts: dict[str, int] = {}
        for fp, label in ((bim, "bim"), (fam, "fam")):
            if fp.exists():
                rows, bad = self._rows_with_columns(fp, 6)
                counts[label] = rows
                if bad:
                    report.invalid.append(f"{bad} row(s) without 6 columns in {fp}")
                else:
                    report.valid += 1

        if bed.exists():
            with open(bed, "rb") as fh:
                magic = fh.read(3)
            if magic != self.PLINK_BED_MAGIC:
                report.invalid.append(f"Invalid .bed magic bytes: {bed}")
            elif "bim" in counts and "fam" in counts:
                expected = 3 + counts["bim"] * math.ceil(counts["fam"] / 4)
                actual = bed.stat().st_size
                if actual != expected:
                    report.invalid.append(
                        f".bed is {actual} bytes; {counts['bim']} variants x "
                        f"{counts['fam']} samples needs {expected}: {bed}"
                    )
                else:
                    report.valid += 1
            else:
                report.valid += 1
        return report

    def validate_batch(self, paths: Sequence[str | Path]) -> ValidationReport:
        """Validate several files.

        ``.vcf`` and ``.vcf.gz`` paths are checked as VCF, anything else as a
        PLINK fileset prefix.
        """
        combined = ValidationReport()
        for p in paths:
            name = str(p)
            is_vcf = name.endswith((".vcf", ".vcf.gz"))
            sub = self.validate_vcf(p) if is_vcf else self.validate_plink_binary(p)
            combined.files_checked += sub.files_checked
            combined.valid += sub.valid
            combined.invalid.extend(sub.invalid)
            combined.warnings.extend(sub.warnings)
        return combined
