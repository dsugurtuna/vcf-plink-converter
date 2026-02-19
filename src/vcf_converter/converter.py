"""Format conversion module.

Wraps PLINK and BCFtools for bidirectional VCF ↔ PLINK conversion
with validation and sample/variant counting.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class ConversionResult:
    """Outcome of a format conversion."""

    input_path: str = ""
    output_prefix: str = ""
    input_format: str = ""
    output_format: str = ""
    sample_count: int = 0
    variant_count: int = 0
    success: bool = False
    message: str = ""


class FormatConverter:
    """Bidirectional VCF ↔ PLINK format converter.

    Parameters
    ----------
    plink_binary : str
        Path to the PLINK 1.9 binary.
    bcftools_binary : str
        Path to the bcftools binary.
    """

    def __init__(
        self,
        plink_binary: str = "plink",
        bcftools_binary: str = "bcftools",
    ) -> None:
        self.plink_binary = plink_binary
        self.bcftools_binary = bcftools_binary

    def _build_vcf_to_plink_cmd(
        self,
        vcf_path: str,
        output_prefix: str,
        extra_args: List[str] | None = None,
    ) -> List[str]:
        cmd = [
            self.plink_binary,
            "--vcf", vcf_path,
            "--make-bed",
            "--out", output_prefix,
            "--allow-extra-chr",
        ]
        if extra_args:
            cmd.extend(extra_args)
        return cmd

    def _build_plink_to_vcf_cmd(
        self,
        bfile_prefix: str,
        output_path: str,
        extra_args: List[str] | None = None,
    ) -> List[str]:
        cmd = [
            self.plink_binary,
            "--bfile", bfile_prefix,
            "--recode", "vcf",
            "--out", output_path,
            "--allow-extra-chr",
        ]
        if extra_args:
            cmd.extend(extra_args)
        return cmd

    def vcf_to_plink(
        self,
        vcf_path: str,
        output_prefix: str,
        extra_args: List[str] | None = None,
    ) -> ConversionResult:
        """Convert VCF to PLINK binary format (.bed/.bim/.fam)."""
        cmd = self._build_vcf_to_plink_cmd(vcf_path, output_prefix, extra_args)
        result = ConversionResult(
            input_path=vcf_path,
            output_prefix=output_prefix,
            input_format="vcf",
            output_format="plink_binary",
        )
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            result.success = True
            result.sample_count = self._count_fam_samples(f"{output_prefix}.fam")
            result.variant_count = self._count_bim_variants(f"{output_prefix}.bim")
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            result.message = str(exc)
        return result

    def plink_to_vcf(
        self,
        bfile_prefix: str,
        output_path: str,
        extra_args: List[str] | None = None,
    ) -> ConversionResult:
        """Convert PLINK binary format to VCF."""
        cmd = self._build_plink_to_vcf_cmd(bfile_prefix, output_path, extra_args)
        result = ConversionResult(
            input_path=bfile_prefix,
            output_prefix=output_path,
            input_format="plink_binary",
            output_format="vcf",
        )
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            result.success = True
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            result.message = str(exc)
        return result

    @staticmethod
    def _count_fam_samples(fam_path: str) -> int:
        """Count samples in a .fam file."""
        path = Path(fam_path)
        if not path.exists():
            return 0
        count = 0
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    count += 1
        return count

    @staticmethod
    def _count_bim_variants(bim_path: str) -> int:
        """Count variants in a .bim file."""
        path = Path(bim_path)
        if not path.exists():
            return 0
        count = 0
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    count += 1
        return count
