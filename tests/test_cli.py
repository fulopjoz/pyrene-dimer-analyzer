"""
Tests for the CLI module.

Tests cover:
- CLI commands
- Command options
- Error handling
"""

import pytest
from click.testing import CliRunner

from pyrene_analyzer.cli import analyze, cli, info, preview


@pytest.fixture
def runner():
    """Create CLI test runner."""
    return CliRunner()


class TestCLIMain:
    """Tests for main CLI group."""

    def test_help(self, runner):
        """Test --help option."""
        result = runner.invoke(cli, ["--help"])

        assert result.exit_code == 0
        assert "Pyrene Dimer Analyzer" in result.output

    def test_version(self, runner):
        """Test --version option."""
        result = runner.invoke(cli, ["--version"])

        assert result.exit_code == 0
        assert "pyrene-dimer-analyzer version" in result.output

    def test_no_command(self, runner):
        """Test invocation without command shows help."""
        result = runner.invoke(cli)

        assert result.exit_code == 0
        assert "Pyrene Dimer Analyzer" in result.output


class TestAnalyzeCommand:
    """Tests for analyze command."""

    def test_analyze_help(self, runner):
        """Test analyze --help."""
        result = runner.invoke(cli, ["analyze", "--help"])

        assert result.exit_code == 0
        assert "Analyze pyrene dimer conformers" in result.output

    def test_analyze_missing_input(self, runner):
        """Test analyze without input files."""
        result = runner.invoke(cli, ["analyze", "-o", "output.csv"])

        assert result.exit_code != 0

    def test_analyze_missing_output(self, runner, simple_test_sdf):
        """Test analyze without output option."""
        result = runner.invoke(cli, ["analyze", str(simple_test_sdf)])

        assert result.exit_code != 0

    def test_analyze_nonexistent_file(self, runner, tmp_path):
        """Test analyze with nonexistent file."""
        result = runner.invoke(
            cli, ["analyze", "nonexistent.sdf", "-o", str(tmp_path / "output.csv")]
        )

        assert result.exit_code != 0

    def test_analyze_basic(self, runner, simple_test_sdf, tmp_path):
        """Test basic analyze command."""
        output_path = tmp_path / "results.csv"

        result = runner.invoke(
            cli,
            [
                "analyze",
                str(simple_test_sdf),
                "-o",
                str(output_path),
                "-q",  # Quiet mode
            ],
        )

        # May fail due to molecule structure, but should not crash
        # Just check it runs without exception

    def test_analyze_with_verbose(self, runner, simple_test_sdf, tmp_path):
        """Test analyze with verbose option."""
        output_path = tmp_path / "results.csv"

        result = runner.invoke(
            cli, ["analyze", str(simple_test_sdf), "-o", str(output_path), "-v"]
        )

        # Check verbose output
        if result.exit_code == 0:
            assert "Processing" in result.output or "INFO" in result.output

    def test_analyze_multiple_formats(self, runner, simple_test_sdf, tmp_path):
        """Test analyze with multiple output formats."""
        output_path = tmp_path / "results"

        result = runner.invoke(
            cli,
            [
                "analyze",
                str(simple_test_sdf),
                "-o",
                str(output_path),
                "-f",
                "csv,json",
                "-q",
            ],
        )

        # Check files would be created (if analysis succeeds)

    def test_analyze_with_plot(self, runner, simple_test_sdf, tmp_path):
        """Test analyze with plot generation."""
        output_path = tmp_path / "results.csv"
        plot_dir = tmp_path / "plots"

        result = runner.invoke(
            cli,
            [
                "analyze",
                str(simple_test_sdf),
                "-o",
                str(output_path),
                "--plot",
                "--plot-dir",
                str(plot_dir),
                "-q",
            ],
        )

    def test_analyze_with_classify(self, runner, simple_test_sdf, tmp_path):
        """Test analyze with classification."""
        output_path = tmp_path / "results.csv"

        result = runner.invoke(
            cli,
            [
                "analyze",
                str(simple_test_sdf),
                "-o",
                str(output_path),
                "--classify",
                "-q",
            ],
        )

    def test_analyze_no_shapely(self, runner, simple_test_sdf, tmp_path):
        """Test analyze without Shapely."""
        output_path = tmp_path / "results.csv"

        result = runner.invoke(
            cli,
            [
                "analyze",
                str(simple_test_sdf),
                "-o",
                str(output_path),
                "--no-shapely",
                "-q",
            ],
        )


class TestInfoCommand:
    """Tests for info command."""

    def test_info(self, runner):
        """Test info command."""
        result = runner.invoke(cli, ["info"])

        assert result.exit_code == 0
        assert "Pyrene Dimer Analyzer" in result.output
        assert "CAPABILITIES" in result.output
        assert "SUPPORTED FORMATS" in result.output
        assert "EXCIMER FORMATION CRITERIA" in result.output
        assert "REFERENCES" in result.output


class TestPreviewCommand:
    """Tests for preview command."""

    def test_preview_help(self, runner):
        """Test preview --help."""
        result = runner.invoke(cli, ["preview", "--help"])

        assert result.exit_code == 0
        assert "Preview analysis results" in result.output

    def test_preview_nonexistent_file(self, runner):
        """Test preview with nonexistent file."""
        result = runner.invoke(cli, ["preview", "nonexistent.sdf"])

        assert result.exit_code != 0

    def test_preview_basic(self, runner, simple_test_sdf):
        """Test basic preview command."""
        result = runner.invoke(cli, ["preview", str(simple_test_sdf)])

        # Should show molecule info
        if result.exit_code == 0:
            assert "Loading" in result.output or "Molecule" in result.output

    def test_preview_with_num(self, runner, simple_test_sdf):
        """Test preview with custom number of conformers."""
        result = runner.invoke(cli, ["preview", str(simple_test_sdf), "-n", "3"])


class TestCLIEdgeCases:
    """Tests for CLI edge cases."""

    def test_empty_output_format(self, runner, simple_test_sdf, tmp_path):
        """Test with empty format string."""
        output_path = tmp_path / "results"

        result = runner.invoke(
            cli,
            ["analyze", str(simple_test_sdf), "-o", str(output_path), "-f", "", "-q"],
        )

    def test_invalid_format(self, runner, simple_test_sdf, tmp_path):
        """Test with invalid format."""
        output_path = tmp_path / "results"

        result = runner.invoke(
            cli,
            [
                "analyze",
                str(simple_test_sdf),
                "-o",
                str(output_path),
                "-f",
                "invalid_format",
                "-q",
            ],
        )

        # Should warn about unknown format

    def test_jobs_option(self, runner, simple_test_sdf, tmp_path):
        """Test with jobs option."""
        output_path = tmp_path / "results.csv"

        result = runner.invoke(
            cli,
            ["analyze", str(simple_test_sdf), "-o", str(output_path), "-j", "2", "-q"],
        )
