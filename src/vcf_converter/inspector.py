"""Lightweight VCF inspection without external tools."""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

_ID = re.compile(r"ID=([^,>]+)")


def iter_vcf_lines(path: str | Path) -> Iterator[str]:
    """Yield lines of a plain or gzip/BGZF-compressed VCF, without newlines."""
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as gz:
            for line in gz:
                yield line.rstrip("\n")
    else:
        with open(path) as fh:
            for line in fh:
                yield line.rstrip("\n")


@dataclass
class InspectionResult:
    """VCF inspection summary."""

    file_path: str = ""
    sample_count: int = 0
    sample_names: list[str] = field(default_factory=list)
    variant_count: int = 0
    contigs: list[str] = field(default_factory=list)
    info_fields: list[str] = field(default_factory=list)
    format_fields: list[str] = field(default_factory=list)
    header_line_count: int = 0


class VCFInspector:
    """Read a VCF header and count data lines."""

    def inspect(self, vcf_path: str | Path) -> InspectionResult:
        result = InspectionResult(file_path=str(vcf_path))
        for line in iter_vcf_lines(vcf_path):
            if line.startswith("##"):
                result.header_line_count += 1
                self._parse_meta_line(line, result)
            elif line.startswith("#CHROM"):
                result.header_line_count += 1
                result.sample_names = line.split("\t")[9:]
                result.sample_count = len(result.sample_names)
            elif line.strip():
                result.variant_count += 1
        return result

    @staticmethod
    def _parse_meta_line(line: str, result: InspectionResult) -> None:
        targets = {
            "##contig=": result.contigs,
            "##INFO=": result.info_fields,
            "##FORMAT=": result.format_fields,
        }
        for prefix, bucket in targets.items():
            if line.startswith(prefix):
                match = _ID.search(line)
                if match:
                    bucket.append(match.group(1))
                return
