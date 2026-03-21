"""Tests for casp16_eval.interface - output parsers (no external tools required)."""

import pytest

from casp16_eval.interface import (
    _parse_ananas_output,
    _parse_ost_interface_output,
    _convert_rna_to_protein_format,
)


class TestParseOstInterfaceOutput:
    SAMPLE = (
        "ICS: 0.72\n"
        "IPS: 0.65\n"
        "i_lDDT: 0.81\n"
        "DockQ: 0.58\n"
    )

    def test_ics_parsed(self):
        result = _parse_ost_interface_output(self.SAMPLE)
        assert result["ics"] == pytest.approx(0.72)

    def test_ips_parsed(self):
        result = _parse_ost_interface_output(self.SAMPLE)
        assert result["ips"] == pytest.approx(0.65)

    def test_i_lddt_parsed(self):
        result = _parse_ost_interface_output(self.SAMPLE)
        assert result["i_lddt"] == pytest.approx(0.81)

    def test_dockq_parsed(self):
        result = _parse_ost_interface_output(self.SAMPLE)
        assert result["dockq"] == pytest.approx(0.58)

    def test_missing_values_return_none(self):
        result = _parse_ost_interface_output("")
        assert result["ics"] is None
        assert result["dockq"] is None


class TestParseAnAnaSOutput:
    SAMPLE = (
        "C1   order=1   RMSD=0.00 Å\n"
        "C2   order=2   RMSD=1.50 Å\n"
        "C3   order=3   RMSD=2.40 Å\n"
        "C4   order=4   RMSD=12.00 Å\n"  # exceeds cutoff=10.0
    )

    def test_highest_order_within_cutoff(self):
        result = _parse_ananas_output(self.SAMPLE, rmsd_cutoff=10.0)
        assert result is not None
        assert result["order"] == 3
        assert result["symmetry"] == "C3"
        assert result["rmsd"] == pytest.approx(2.40)

    def test_all_exceed_cutoff_returns_none(self):
        output = "C2   order=2   RMSD=15.0 Å\n"
        result = _parse_ananas_output(output, rmsd_cutoff=10.0)
        assert result is None

    def test_empty_output_returns_none(self):
        result = _parse_ananas_output("", rmsd_cutoff=10.0)
        assert result is None


class TestConvertRnaToProteinFormat:
    def test_phosphorus_renamed_to_ca(self, tmp_path):
        pdb_content = (
            "ATOM      1  P     G A   1       1.000   2.000   3.000  1.00  0.00           P\n"
        )
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(pdb_content)

        out = _convert_rna_to_protein_format(str(pdb_file), str(tmp_path))
        assert out is not None
        content = open(out).read()
        assert " CA " in content

    def test_residue_renamed_to_ala(self, tmp_path):
        pdb_content = (
            "ATOM      1  P     G A   1       1.000   2.000   3.000  1.00  0.00           P\n"
        )
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(pdb_content)

        out = _convert_rna_to_protein_format(str(pdb_file), str(tmp_path))
        assert out is not None
        content = open(out).read()
        assert "ALA" in content

    def test_non_rna_atom_unchanged(self, tmp_path):
        pdb_content = (
            "ATOM      2  CA  ALA A   2       4.000   5.000   6.000  1.00  0.00           C\n"
        )
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(pdb_content)

        out = _convert_rna_to_protein_format(str(pdb_file), str(tmp_path))
        assert out is not None
        content = open(out).read()
        # ALA residue should be preserved unchanged
        assert " CA  ALA" in content

    def test_nonexistent_file_returns_none(self, tmp_path):
        result = _convert_rna_to_protein_format("/nonexistent/path.pdb", str(tmp_path))
        assert result is None
