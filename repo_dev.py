#!/usr/bin/python
"""
Simple shim for running the repo program during development.
"""
import ga4gh.cli
# PYTHON_ARGCOMPLETE_OK

if __name__ == "__main__":
    ga4gh.cli.repo_main()
