"""
Setup script for pyrene-dimer-analyzer package.

Install with:
    pip install .
    pip install -e .  # Development mode
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_path = Path(__file__).parent / "README.md"
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")
else:
    long_description = "Automated geometric analysis of pyrene dimer conformers"

# Read version from package
version = {}
version_path = Path(__file__).parent / "pyrene_analyzer" / "__init__.py"
if version_path.exists():
    with open(version_path) as f:
        for line in f:
            if line.startswith("__version__"):
                exec(line, version)
                break

setup(
    name="pyrene-dimer-analyzer",
    version=version.get("__version__", "1.0.0"),
    description="Automated geometric analysis of pyrene dimer conformers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Research Team, UCT Prague",
    author_email="research@vscht.cz",
    url="https://github.com/research-team/pyrene-dimer-analyzer",
    project_urls={
        "Documentation": "https://pyrene-dimer-analyzer.readthedocs.io/",
        "Source": "https://github.com/research-team/pyrene-dimer-analyzer",
        "Bug Tracker": "https://github.com/research-team/pyrene-dimer-analyzer/issues",
    },
    packages=find_packages(exclude=["tests", "tests.*", "docs", "examples", "notebooks"]),
    python_requires=">=3.9",
    install_requires=[
        "rdkit>=2023.9.1",
        "numpy>=1.24.0",
        "scipy>=1.11.0",
        "pandas>=2.0.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "shapely>=2.0.0",
        "click>=8.1.0",
        "tqdm>=4.65.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.4.0",
            "isort>=5.12.0",
        ],
        "docs": [
            "sphinx>=7.0.0",
            "sphinx-rtd-theme>=1.3.0",
            "myst-parser>=2.0.0",
        ],
        "notebooks": [
            "jupyter>=1.0.0",
            "nbconvert>=7.0.0",
            "ipywidgets>=8.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "pyrene-analyze=pyrene_analyzer.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords=[
        "pyrene",
        "excimer",
        "conformational analysis",
        "molecular geometry",
        "rdkit",
        "cheminformatics",
        "computational chemistry",
    ],
    include_package_data=True,
    zip_safe=False,
)
