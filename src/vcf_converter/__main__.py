"""Command-line interface: ``python -m vcf_converter <command>``."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from .converter import FormatConverter
from .inspector import VCFInspector
from .validator import FileValidator


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="vcf_converter")
    sub = parser.add_subparsers(dest="command", required=True)
    inspect = sub.add_parser("inspect", help="summarise a VCF header and counts")
    inspect.add_argument("vcf")
    validate = sub.add_parser("validate", help="check VCFs or PLINK prefixes")
    validate.add_argument("paths", nargs="+")
    to_plink = sub.add_parser("vcf-to-plink", help="VCF -> .bed/.bim/.fam (PLINK)")
    to_plink.add_argument("vcf")
    to_plink.add_argument("out_prefix")
    to_vcf = sub.add_parser("plink-to-vcf", help=".bed/.bim/.fam -> VCF (PLINK)")
    to_vcf.add_argument("bfile")
    to_vcf.add_argument("out", help="prefix, name.vcf or name.vcf.gz")
    for p in (to_plink, to_vcf):
        p.add_argument("--plink", default="plink", help="PLINK 1.9 executable")
    args = parser.parse_args(argv)

    if args.command == "inspect":
        info = VCFInspector().inspect(args.vcf)
        print(f"samples:  {info.sample_count} ({', '.join(info.sample_names[:5])})")
        print(f"variants: {info.variant_count}")
        print(f"contigs:  {', '.join(info.contigs) or '-'}")
        print(f"INFO:     {', '.join(info.info_fields) or '-'}")
        print(f"FORMAT:   {', '.join(info.format_fields) or '-'}")
        return 0
    if args.command == "validate":
        report = FileValidator().validate_batch(args.paths)
        for problem in report.invalid:
            print(f"INVALID  {problem}")
        print(f"{report.valid} checks passed, {len(report.invalid)} problems")
        return 0 if report.all_valid else 1

    converter = FormatConverter(plink_binary=args.plink)
    if args.command == "vcf-to-plink":
        result = converter.vcf_to_plink(args.vcf, args.out_prefix)
    else:
        result = converter.plink_to_vcf(args.bfile, args.out)
    if not result.success:
        print(f"failed: {result.message}")
        return 1
    print(
        f"wrote {result.output_path}: {result.sample_count} samples, "
        f"{result.variant_count} variants"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
