#!/usr/bin/env python3
"""
Automated Workflow Example wrapper.

This file simply forwards to the production-ready workflow in
`automated_workflow.py`. Keep using this script for demonstrations or
quick smoke tests, but the canonical entrypoint is now the non-example
module.
"""

from automated_workflow import main


if __name__ == "__main__":
    main()
