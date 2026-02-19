"""VCF-PLINK Converter — bidirectional VCF and PLINK format conversion."""

__version__ = "1.0.0"

from .converter import FormatConverter, ConversionResult
from .validator import FileValidator, ValidationReport
from .inspector import VCFInspector, InspectionResult

__all__ = [
    "FormatConverter",
    "ConversionResult",
    "FileValidator",
    "ValidationReport",
    "VCFInspector",
    "InspectionResult",
]
