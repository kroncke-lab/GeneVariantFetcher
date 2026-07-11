"""Acquisition expected-value (EV) package.

Exposes the abstract-only, gold-free, pure-stdlib yield scorer in
``predict_yield`` as an importable module so the source-acquisition worklist can
rank the un-downloaded fetch tail by predicted per-paper carrier/variant yield.
Importing this package (and ``predict_yield``) has no side effects and makes no
network or LLM calls.
"""
