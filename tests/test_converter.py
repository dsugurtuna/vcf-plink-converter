"""Tests for FormatConverter and FileValidator."""

from vcf_converter.converter import FormatConverter, ConversionResult
from vcf_converter.validator import FileValidator


class TestFormatConverter:
    def test_build_vcf_to_plink_cmd(self):
        conv = FormatConverter()
        cmd = conv._build_vcf_to_plink_cmd("input.vcf", "output")
        assert "--vcf" in cmd
        assert "input.vcf" in cmd
        assert "--make-bed" in cmd

    def test_build_plink_to_vcf_cmd(self):
        conv = FormatConverter()
        cmd = conv._build_plink_to_vcf_cmd("input", "output")
        assert "--bfile" in cmd
        assert "--recode" in cmd
        assert "vcf" in cmd

    def test_extra_args(self):
        conv = FormatConverter()
        cmd = conv._build_vcf_to_plink_cmd("in.vcf", "out", ["--maf", "0.01"])
        assert "--maf" in cmd
        assert "0.01" in cmd

    def test_count_fam_samples(self, tmp_path):
        fam = tmp_path / "test.fam"
        fam.write_text("FAM1 IND1 0 0 1 -9\nFAM2 IND2 0 0 2 -9\n")
        count = FormatConverter._count_fam_samples(str(fam))
        assert count == 2

    def test_count_bim_variants(self, tmp_path):
        bim = tmp_path / "test.bim"
        bim.write_text("1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tC\tT\n1\trs3\t0\t300\tG\tA\n")
        count = FormatConverter._count_bim_variants(str(bim))
        assert count == 3

    def test_count_missing_fam(self):
        assert FormatConverter._count_fam_samples("/nonexistent.fam") == 0

    def test_conversion_result_defaults(self):
        result = ConversionResult()
        assert not result.success
        assert result.sample_count == 0


class TestFileValidator:
    def test_validate_vcf_valid(self, tmp_path):
        vcf = tmp_path / "test.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\n1\t100\trs1\n")
        val = FileValidator()
        report = val.validate_vcf(vcf)
        assert report.all_valid

    def test_validate_vcf_invalid(self, tmp_path):
        vcf = tmp_path / "test.vcf"
        vcf.write_text("not a vcf file\n")
        val = FileValidator()
        report = val.validate_vcf(vcf)
        assert not report.all_valid

    def test_validate_vcf_missing(self):
        val = FileValidator()
        report = val.validate_vcf("/nonexistent.vcf")
        assert not report.all_valid

    def test_validate_plink_binary(self, tmp_path):
        prefix = tmp_path / "test"
        (tmp_path / "test.bed").write_bytes(b"\x6c\x1b\x01" + b"\x00" * 10)
        (tmp_path / "test.bim").write_text("1\trs1\t0\t100\tA\tG\n")
        (tmp_path / "test.fam").write_text("FAM1 IND1 0 0 1 -9\n")
        val = FileValidator()
        report = val.validate_plink_binary(str(prefix))
        assert report.all_valid
        assert report.valid == 3

    def test_validate_plink_missing_bed(self, tmp_path):
        prefix = tmp_path / "test"
        (tmp_path / "test.bim").write_text("data\n")
        (tmp_path / "test.fam").write_text("data\n")
        val = FileValidator()
        report = val.validate_plink_binary(str(prefix))
        assert not report.all_valid
