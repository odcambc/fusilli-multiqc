"""Tests for fusilli_multiqc.utils."""

import tempfile
from pathlib import Path

import pandas as pd

from fusilli_multiqc.utils import parse_csv_file, parse_json_file


class TestParseCsvFile:
    def test_valid_csv(self, tmp_path):
        csv = tmp_path / "data.csv"
        csv.write_text("a,b\n1,2\n3,4\n")
        df = parse_csv_file(str(csv))
        assert df is not None
        assert list(df.columns) == ["a", "b"]
        assert len(df) == 2

    def test_missing_file(self):
        assert parse_csv_file("/nonexistent/file.csv") is None

    def test_empty_csv(self, tmp_path):
        csv = tmp_path / "empty.csv"
        csv.write_text("a,b\n")
        assert parse_csv_file(str(csv)) is None

    def test_malformed_csv(self, tmp_path):
        csv = tmp_path / "bad.csv"
        csv.write_text("this is not csv\x00\x01\x02")
        # Should not raise; returns a DataFrame or None
        result = parse_csv_file(str(csv))
        # Malformed but parseable as single-column — either result is fine
        assert result is None or isinstance(result, pd.DataFrame)


class TestParseJsonFile:
    def test_valid_json(self, tmp_path):
        jf = tmp_path / "data.json"
        jf.write_text('{"key": "value"}')
        result = parse_json_file(str(jf))
        assert result == {"key": "value"}

    def test_missing_json(self):
        assert parse_json_file("/nonexistent/file.json") is None

    def test_invalid_json(self, tmp_path):
        jf = tmp_path / "bad.json"
        jf.write_text("not json at all")
        assert parse_json_file(str(jf)) is None
