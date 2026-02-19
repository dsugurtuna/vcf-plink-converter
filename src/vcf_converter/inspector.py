"""VCF inspection module.

Provides lightweight VCF header and content inspection without
requiring external tools.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


@dataclass
class InspectionResult:
    """VCF inspection summary."""

    file_path: str = ""
    sample_count: int = 0
    variant_count: int = 0
    contigs: List[str] = field(default_factory=list)
    info_fields: List[str] = field(default_factory=list)
    format_fields: List[str] = field(default_factory=list)
    header_line_count: int = 0


class VCFInspector:
    """Inspect VCF files to extract metadata and counts.

    Parses the VCF header to extract sample names, contig info,
    and INFO/FORMAT fields, then counts data lines for variant total.
    """

    def inspect(self, vcf_path: str | Path) -> InspectionResult:
        """Inspect a VCF file."""
        path = Path(vcf_path)
        result = InspectionResult(file_path=str(vcf_path))

        opener = gzip.open if path.suffix == ".gz" else open
        mode = "rt" if path.suffix == ".gz" else "r"

        with opener(path, mode) as fh:  # type: ignore[call-overload]
            for line in fh:
                line = line.strip()
                if line.startswith("##"):
                    result.header_line_count += 1
                    self._parse_meta_line(line, result)
                elif line.startswith("#CHROM"):
                    result.header_line_count += 1
                    cols = line.split("\t")
                    if len(cols) > 9:
                        result.sample_count = len(cols) - 9
                else:
                    if line:
                        result.variant_count += 1

        return result

    @staticmethod
    def _parse_meta_line(line: str, result: InspectionResult) -> None:
        if line.startswith("##contig="):
            # Extract contig ID
            start = line.find("ID=")
            if start >= 0:
                end = line.find(",", start)
                if end < 0:
                    end = line.find(">", start)
                contig_id = line[start + 3: end]
                if contig_id:
                    result.contigs.append(contig_id)
        elif line.startswith("##INFO="):
            start = line.find("ID=")
            if start >= 0:
                end = line.find(",", start)
                field_id = line[start + 3: end]
                if field_id:
                    result.info_fields.append(field_id)
        elif line.startswith("##FORMAT="):
            start = line.find("ID=")
            if start >= 0:
                end = line.find(",", start)
                field_id = line[start + 3: end]
                if field_id:
                    result.format_fields.append(field_id)
