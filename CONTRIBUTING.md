# Contributing to WGS Analysis Pipeline

Thank you for your interest in contributing to the WGS Analysis Pipeline! We welcome contributions from the community.

## How to Contribute

### Reporting Issues

- Use the GitHub issue tracker to report bugs or request features
- Provide detailed information including:
  - Steps to reproduce the issue
  - Expected vs. actual behavior
  - Your environment (OS, conda version, etc.)
  - Relevant log files or error messages

### Contributing Code

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Make your changes
4. Run tests and linting: `snakemake --lint` for Snakemake files
5. Commit your changes: `git commit -m "Add your message"`
6. Push to your fork: `git push origin feature/your-feature-name`
7. Create a Pull Request

### Code Style

- Follow existing code conventions
- Use descriptive variable and function names
- Add comments for complex logic
- Keep lines under 100 characters where possible

### Testing

- Test your changes with sample data
- Ensure Snakemake workflows run without errors
- Verify output files are generated correctly

## Development Setup

1. Clone the repository
2. Create the conda environment: `conda env create -f environment.yml`
3. Activate: `conda activate wgs_analysis`
4. Make your changes
5. Test thoroughly

## Documentation

- Update README.md for new features
- Add docstrings to new functions
- Update workflow documentation in `docs/`

## Questions?

If you have questions about contributing, please open an issue or contact the maintainers.