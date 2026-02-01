"""
Tests for the visualization module.

Tests cover:
- Angle vs energy plots
- Distance vs overlap plots
- Conformer distribution plots
- Variant comparison plots
- Summary figures
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for testing

import matplotlib.pyplot as plt
import pytest

from pyrene_analyzer.visualization import (
    plot_angle_vs_energy,
    plot_distance_vs_overlap,
    plot_conformer_distribution,
    plot_variant_comparison,
    plot_energy_landscape,
    plot_correlation_matrix,
    create_summary_figure,
)


class TestPlotAngleVsEnergy:
    """Tests for plot_angle_vs_energy function."""
    
    def test_basic_plot(self, sample_analysis_df):
        """Test basic angle vs energy plot."""
        fig = plot_angle_vs_energy(sample_analysis_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, sample_analysis_df, tmp_path):
        """Test saving plot to file."""
        output_path = tmp_path / "angle_energy.png"
        fig = plot_angle_vs_energy(sample_analysis_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)
    
    def test_plot_with_color_by(self, multi_variant_df):
        """Test plot with color by molecule."""
        fig = plot_angle_vs_energy(multi_variant_df, color_by='molecule')
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_with_custom_title(self, sample_analysis_df):
        """Test plot with custom title."""
        fig = plot_angle_vs_energy(
            sample_analysis_df, 
            title="Custom Title"
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_without_excimer_threshold(self, sample_analysis_df):
        """Test plot without excimer threshold lines."""
        fig = plot_angle_vs_energy(
            sample_analysis_df,
            show_excimer_threshold=False
        )
        
        assert fig is not None
        plt.close(fig)


class TestPlotDistanceVsOverlap:
    """Tests for plot_distance_vs_overlap function."""
    
    def test_basic_plot(self, sample_analysis_df):
        """Test basic distance vs overlap plot."""
        fig = plot_distance_vs_overlap(sample_analysis_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, sample_analysis_df, tmp_path):
        """Test saving plot to file."""
        output_path = tmp_path / "distance_overlap.png"
        fig = plot_distance_vs_overlap(sample_analysis_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)
    
    def test_plot_with_color_by(self, multi_variant_df):
        """Test plot with color by molecule."""
        fig = plot_distance_vs_overlap(multi_variant_df, color_by='molecule')
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_without_optimal_region(self, sample_analysis_df):
        """Test plot without optimal region highlight."""
        fig = plot_distance_vs_overlap(
            sample_analysis_df,
            show_optimal_region=False
        )
        
        assert fig is not None
        plt.close(fig)


class TestPlotConformerDistribution:
    """Tests for plot_conformer_distribution function."""
    
    def test_basic_plot(self, sample_analysis_df):
        """Test basic distribution plot."""
        fig = plot_conformer_distribution(sample_analysis_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, sample_analysis_df, tmp_path):
        """Test saving plot to file."""
        output_path = tmp_path / "distributions.png"
        fig = plot_conformer_distribution(sample_analysis_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)
    
    def test_plot_custom_parameters(self, sample_analysis_df):
        """Test plot with custom parameters."""
        fig = plot_conformer_distribution(
            sample_analysis_df,
            parameters=['plane_angle_deg', 'interplane_distance_A']
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_single_parameter(self, sample_analysis_df):
        """Test plot with single parameter."""
        fig = plot_conformer_distribution(
            sample_analysis_df,
            parameters=['plane_angle_deg']
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_invalid_parameters(self, sample_analysis_df):
        """Test plot with invalid parameters."""
        with pytest.raises(ValueError):
            plot_conformer_distribution(
                sample_analysis_df,
                parameters=['nonexistent_column']
            )


class TestPlotVariantComparison:
    """Tests for plot_variant_comparison function."""
    
    def test_basic_plot(self, multi_variant_df):
        """Test basic variant comparison plot."""
        fig = plot_variant_comparison(multi_variant_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, multi_variant_df, tmp_path):
        """Test saving plot to file."""
        output_path = tmp_path / "comparison.png"
        fig = plot_variant_comparison(multi_variant_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)
    
    def test_plot_specific_variants(self, multi_variant_df):
        """Test plot with specific variants."""
        fig = plot_variant_comparison(
            multi_variant_df,
            variants=['Et', 'iPr']
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_different_parameter(self, multi_variant_df):
        """Test plot with different parameter."""
        fig = plot_variant_comparison(
            multi_variant_df,
            parameter='interplane_distance_A'
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_plot_invalid_variants(self, multi_variant_df):
        """Test plot with invalid variants."""
        with pytest.raises(ValueError):
            plot_variant_comparison(
                multi_variant_df,
                variants=['nonexistent_variant']
            )


class TestPlotEnergyLandscape:
    """Tests for plot_energy_landscape function."""
    
    def test_basic_plot(self, sample_analysis_df):
        """Test basic energy landscape plot."""
        fig = plot_energy_landscape(sample_analysis_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, sample_analysis_df, tmp_path):
        """Test saving plot to file."""
        output_path = tmp_path / "landscape.png"
        fig = plot_energy_landscape(sample_analysis_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)
    
    def test_plot_no_energy(self, multi_variant_df):
        """Test plot with no energy data."""
        with pytest.raises(ValueError):
            plot_energy_landscape(multi_variant_df)


class TestPlotCorrelationMatrix:
    """Tests for plot_correlation_matrix function."""
    
    def test_basic_plot(self, sample_analysis_df):
        """Test basic correlation matrix plot."""
        fig = plot_correlation_matrix(sample_analysis_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, sample_analysis_df, tmp_path):
        """Test saving plot to file."""
        output_path = tmp_path / "correlation.png"
        fig = plot_correlation_matrix(sample_analysis_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)


class TestCreateSummaryFigure:
    """Tests for create_summary_figure function."""
    
    def test_basic_plot(self, sample_analysis_df):
        """Test basic summary figure."""
        fig = create_summary_figure(sample_analysis_df)
        
        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_plot_with_output(self, sample_analysis_df, tmp_path):
        """Test saving summary figure to file."""
        output_path = tmp_path / "summary.png"
        fig = create_summary_figure(sample_analysis_df, output=output_path)
        
        assert output_path.exists()
        plt.close(fig)
    
    def test_plot_with_classification(self, sample_analysis_df):
        """Test summary figure with classification column."""
        from pyrene_analyzer.core import PyreneDimerAnalyzer
        
        analyzer = PyreneDimerAnalyzer(verbose=False)
        df = analyzer.add_classification(sample_analysis_df)
        
        fig = create_summary_figure(df)
        
        assert fig is not None
        plt.close(fig)


class TestPlotStyling:
    """Tests for plot styling and formatting."""
    
    def test_custom_figsize(self, sample_analysis_df):
        """Test custom figure size."""
        fig = plot_angle_vs_energy(
            sample_analysis_df,
            figsize=(12, 8)
        )
        
        assert fig.get_figwidth() == 12
        assert fig.get_figheight() == 8
        plt.close(fig)
    
    def test_custom_alpha(self, sample_analysis_df):
        """Test custom alpha transparency."""
        fig = plot_angle_vs_energy(
            sample_analysis_df,
            alpha=0.5
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_custom_marker_size(self, sample_analysis_df):
        """Test custom marker size."""
        fig = plot_angle_vs_energy(
            sample_analysis_df,
            marker_size=100
        )
        
        assert fig is not None
        plt.close(fig)
