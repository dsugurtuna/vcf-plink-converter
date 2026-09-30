"""VCF-PLINK Converter — bidirectional VCF and PLINK format conversion."""

__version__ = "1.0.0"

from .converter import ConversionResult, FormatConverter
from .inspector import InspectionResult, VCFInspector
from .validator import FileValidator, ValidationReport

__all__ = [
    "FormatConverter",
    "ConversionResult",
    "FileValidator",
    "ValidationReport",
    "VCFInspector",
    "InspectionResult",
]
