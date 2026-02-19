"""Tests for VCFInspector."""

from vcf_converter.inspector import VCFInspector


class TestVCFInspector:
    def _write_vcf(self, tmp_path, content):
        vcf = tmp_path / "test.vcf"
        vcf.write_text(content)
        return vcf

    def test_basic_inspection(self, tmp_path):
        content = (
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "##contig=<ID=1,length=249250621>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            "1\t100\trs1\tA\tG\t.\tPASS\tAC=1\tGT\t0/1\n"
            "1\t200\trs2\tC\tT\t.\tPASS\tAC=2\tGT\t1/1\n"
        )
        vcf = self._write_vcf(tmp_path, content)
        inspector = VCFInspector()
        result = inspector.inspect(vcf)
        assert result.sample_count == 1
        assert result.variant_count == 2
        assert "1" in result.contigs
        assert "AC" in result.info_fields
        assert "GT" in result.format_fields

    def test_multiple_samples(self, tmp_path):
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n"
        )
        vcf = self._write_vcf(tmp_path, content)
        result = VCFInspector().inspect(vcf)
        assert result.sample_count == 3
        assert result.variant_count == 1

    def test_empty_vcf(self, tmp_path):
        content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        vcf = self._write_vcf(tmp_path, content)
        result = VCFInspector().inspect(vcf)
        assert result.sample_count == 0
        assert result.variant_count == 0

    def test_header_line_count(self, tmp_path):
        content = (
            "##fileformat=VCFv4.2\n"
            "##source=test\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        vcf = self._write_vcf(tmp_path, content)
        result = VCFInspector().inspect(vcf)
        assert result.header_line_count == 3
