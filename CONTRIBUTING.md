# Contributing to pyrene-dimer-analyzer

Thank you for your interest in contributing to `pyrene-dimer-analyzer`! This document provides guidelines and instructions for contributing.

## Code of Conduct

This project adheres to a Code of Conduct. By participating, you are expected to uphold this code.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with:

1. A clear, descriptive title
2. Steps to reproduce the issue
3. Expected behavior
4. Actual behavior
5. Your environment (Python version, OS, package versions)

### Suggesting Features

Feature requests are welcome! Please open an issue with:

1. A clear description of the feature
2. The use case or problem it solves
3. Any relevant examples or references

### Pull Requests

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature`)
3. Make your changes
4. Run tests (`pytest tests/`)
5. Format code (`black pyrene_analyzer tests`)
6. Commit your changes (`git commit -m 'Add feature'`)
7. Push to the branch (`git push origin feature/your-feature`)
8. Open a Pull Request

## Development Setup

```bash
# Clone your fork
git clone https://github.com/your-username/pyrene-dimer-analyzer
cd pyrene-dimer-analyzer

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest tests/ -v --cov=pyrene_analyzer
```

## Code Style

We use the following tools to maintain code quality:

- **Black** for code formatting
- **isort** for import sorting
- **flake8** for linting
- **mypy** for type checking

Before submitting a PR, please run:

```bash
black pyrene_analyzer tests
isort pyrene_analyzer tests
flake8 pyrene_analyzer
```

## Testing

All new features should include tests. We aim for >80% code coverage.

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=pyrene_analyzer --cov-report=html
```

## Documentation

Please update documentation for any new features or changes:

1. Update docstrings in the code
2. Update relevant `.md` files in `docs/`
3. Add examples if appropriate

## Questions?

Feel free to open an issue for any questions about contributing.
