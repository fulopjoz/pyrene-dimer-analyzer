# Installation

This guide provides instructions for installing `pyrene-dimer-analyzer` and its dependencies.

## Recommended Installation (from PyPI)

The easiest way to install the package is from the Python Package Index (PyPI) using `pip`:

```bash
pip install pyrene-dimer-analyzer
```

This will download and install the latest stable version of the package along with all required dependencies.

## Installation from Source

If you want to install the latest development version, you can clone the repository from GitHub and install it in editable mode:

```bash
# 1. Clone the repository
git clone https://github.com/research-team/pyrene-dimer-analyzer

# 2. Navigate to the project directory
cd pyrene-dimer-analyzer

# 3. Install in editable mode
pip install -e .
```

Installing in editable mode (`-e`) means that any changes you make to the source code will be immediately available without needing to reinstall the package.

## Dependencies

`pyrene-dimer-analyzer` relies on several powerful open-source libraries:

-   **RDKit**: For core cheminformatics functionality.
-   **NumPy & SciPy**: For numerical and scientific computing.
-   **Pandas**: For data manipulation and analysis.
-   **Matplotlib & Seaborn**: For data visualization.
-   **Shapely**: For accurate π-overlap calculations.
-   **Click**: For the command-line interface.
-   **tqdm**: For progress bars.

All of these dependencies will be automatically installed when you install `pyrene-dimer-analyzer` using `pip`.

## Development Installation

If you plan to contribute to the development of the package, you should install the development dependencies as well:

```bash
# After cloning the repository
pip install -e ".[dev,docs,notebooks]"
```

This will install all the tools required for testing, linting, formatting, and building the documentation.

## Verifying the Installation

After installation, you can verify that the package is working correctly by running the following command:

```bash
pyrene-analyze --version
```

This should print the installed version of the package. You can also get more detailed information by running:

```bash
pyrene-analyze info
```
