"""File validation module.

Validates VCF and PLINK binary file integrity before conversion.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import List


@dataclass
class ValidationReport:
    """Validation outcome for a set of files."""

    files_checked: int = 0
    valid: int = 0
    invalid: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def all_valid(self) -> bool:
        return len(self.invalid) == 0


class FileValidator:
    """Validate VCF and PLINK binary files.

    Performs lightweight structural checks — header presence,
    magic bytes, and file triads.
    """

    VCF_HEADER = "##fileformat=VCF"
    PLINK_BED_MAGIC = b"\x6c\x1b\x01"

    def validate_vcf(self, vcf_path: str | Path) -> ValidationReport:
        """Validate a VCF file by checking the header line."""
        report = ValidationReport(files_checked=1)
        path = Path(vcf_path)

        if not path.exists():
            report.invalid.append(f"File not found: {vcf_path}")
            return report

        try:
            if path.suffix == ".gz":
                with gzip.open(path, "rt") as fh:
                    first_line = fh.readline().strip()
            else:
                with open(path) as fh:
                    first_line = fh.readline().strip()

            if first_line.startswith(self.VCF_HEADER):
                report.valid = 1
            else:
                report.invalid.append(f"Missing VCF header: {vcf_path}")
        except (OSError, gzip.BadGzipFile) as exc:
            report.invalid.append(f"Read error: {exc}")

        return report

    def validate_plink_binary(self, bfile_prefix: str | Path) -> ValidationReport:
        """Validate a PLINK binary fileset (.bed/.bim/.fam)."""
        prefix = Path(bfile_prefix)
        bed = Path(f"{prefix}.bed")
        bim = Path(f"{prefix}.bim")
        fam = Path(f"{prefix}.fam")
        report = ValidationReport(files_checked=3)

        for fp, label in [(bed, "bed"), (bim, "bim"), (fam, "fam")]:
            if not fp.exists():
                report.invalid.append(f"Missing .{label} file: {fp}")

        if bed.exists():
            with open(bed, "rb") as fh:
                magic = fh.read(3)
            if magic == self.PLINK_BED_MAGIC:
                report.valid += 1
            else:
                report.invalid.append(f"Invalid .bed magic bytes: {bed}")
        if bim.exists():
            report.valid += 1
        if fam.exists():
            report.valid += 1

        return report

    def validate_batch(self, paths: List[str | Path]) -> ValidationReport:
        """Validate multiple files (VCF or PLINK prefix detection)."""
        combined = ValidationReport()
        for p in paths:
            path = Path(p)
            if path.suffix in (".vcf", ".gz"):
                sub = self.validate_vcf(p)
            else:
                sub = self.validate_plink_binary(p)
            combined.files_checked += sub.files_checked
            combined.valid += sub.valid
            combined.invalid.extend(sub.invalid)
            combined.warnings.extend(sub.warnings)
        return combined
