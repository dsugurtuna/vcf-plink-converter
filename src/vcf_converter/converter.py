"""VCF <-> PLINK binary conversion with PLINK 1.9.

PLINK 1.9's defaults lose information in a round trip, so they are
overridden:

* ``--double-id`` on import: otherwise a VCF sample ID with one underscore
  is split into FID and IID (``SAMPLE_001`` -> ``SAMPLE`` / ``001``), and
  an ID with two underscores is rejected.
* ``--recode vcf-iid`` on export: otherwise sample names are written as
  ``FID_IID``.
* ``--keep-allele-order`` both ways: otherwise PLINK may make the major
  allele A2, and A2 is written as REF on export.
"""

from __future__ import annotations

import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .inspector import VCFInspector

Runner = Callable[[list[str]], "subprocess.CompletedProcess[str]"]


def subprocess_runner(argv: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(argv, capture_output=True, text=True, check=False)


@dataclass
class ConversionResult:
    """Outcome of a format conversion."""

    input_path: str = ""
    output_prefix: str = ""
    input_format: str = ""
    output_format: str = ""
    output_path: str = ""
    sample_count: int = 0
    variant_count: int = 0
    success: bool = False
    message: str = ""


class FormatConverter:
    """Convert between VCF and PLINK binary filesets.

    Parameters
    ----------
    plink_binary : str
        PLINK 1.9 executable.
    runner : callable, optional
        Runs one command and returns a CompletedProcess (for testing).
    """

    def __init__(
        self, plink_binary: str = "plink", runner: Runner | None = None
    ) -> None:
        self.plink_binary = plink_binary
        self.runner = runner or subprocess_runner

    def _build_vcf_to_plink_cmd(
        self, vcf_path: str, output_prefix: str, extra_args: list[str] | None = None
    ) -> list[str]:
        return [
            self.plink_binary,
            "--vcf",
            vcf_path,
            "--double-id",
            "--keep-allele-order",
            "--allow-extra-chr",
            "--make-bed",
            "--out",
            output_prefix,
            *(extra_args or []),
        ]

    def _build_plink_to_vcf_cmd(
        self,
        bfile_prefix: str,
        output_prefix: str,
        extra_args: list[str] | None = None,
        compress: bool = False,
    ) -> list[str]:
        return [
            self.plink_binary,
            "--bfile",
            bfile_prefix,
            "--recode",
            "vcf-iid",
            *(["bgz"] if compress else []),
            "--keep-allele-order",
            "--allow-extra-chr",
            "--out",
            output_prefix,
            *(extra_args or []),
        ]

    def _run(self, cmd: list[str], result: ConversionResult) -> bool:
        try:
            proc = self.runner(cmd)
        except FileNotFoundError as exc:  # PLINK not installed
            result.message = str(exc)
            return False
        if proc.returncode != 0:
            result.message = (proc.stderr or proc.stdout).strip()
            return False
        return True

    def vcf_to_plink(
        self, vcf_path: str, output_prefix: str, extra_args: list[str] | None = None
    ) -> ConversionResult:
        """Write ``<output_prefix>.bed/.bim/.fam`` from a VCF (plain or gzipped)."""
        result = ConversionResult(
            input_path=vcf_path,
            output_prefix=output_prefix,
            input_format="vcf",
            output_format="plink_binary",
            output_path=output_prefix,
        )
        cmd = self._build_vcf_to_plink_cmd(vcf_path, output_prefix, extra_args)
        if self._run(cmd, result):
            result.success = True
            result.sample_count = self._count_fam_samples(f"{output_prefix}.fam")
            result.variant_count = self._count_bim_variants(f"{output_prefix}.bim")
        return result

    def plink_to_vcf(
        self, bfile_prefix: str, output_path: str, extra_args: list[str] | None = None
    ) -> ConversionResult:
        """Write a VCF from a PLINK fileset.

        ``output_path`` may be a prefix, ``name.vcf`` or ``name.vcf.gz``
        (block-gzipped). PLINK adds the extension itself.
        """
        prefix, compress = output_path, False
        for suffix, gz in ((".vcf.gz", True), (".vcf", False)):
            if output_path.endswith(suffix):
                prefix, compress = output_path[: -len(suffix)], gz
                break
        written = prefix + (".vcf.gz" if compress else ".vcf")
        result = ConversionResult(
            input_path=bfile_prefix,
            output_prefix=prefix,
            input_format="plink_binary",
            output_format="vcf",
            output_path=written,
        )
        cmd = self._build_plink_to_vcf_cmd(bfile_prefix, prefix, extra_args, compress)
        if self._run(cmd, result):
            result.success = True
            if Path(written).exists():
                info = VCFInspector().inspect(written)
                result.sample_count = info.sample_count
                result.variant_count = info.variant_count
        return result

    @staticmethod
    def _count_lines(path: str) -> int:
        p = Path(path)
        if not p.exists():
            return 0
        with open(p) as fh:
            return sum(1 for line in fh if line.strip())

    @staticmethod
    def _count_fam_samples(fam_path: str) -> int:
        """Samples in a .fam file."""
        return FormatConverter._count_lines(fam_path)

    @staticmethod
    def _count_bim_variants(bim_path: str) -> int:
        """Variants in a .bim file."""
        return FormatConverter._count_lines(bim_path)
