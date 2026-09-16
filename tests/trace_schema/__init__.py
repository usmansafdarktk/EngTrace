"""Structured trace nodes (Phase 3, D3.3).

The specification lives in `docs/re-implementation-sep/track-a/phase3_node_types.md`.
This package holds only the tooling: an extractor that lifts a template's
`trace_nodes` local out of a generation call, and a corpus dumper that writes
node/prose pairs for a verifier to be written against.

Nothing here re-implements a template. The node a template builds IS the object
its prose is rendered from, so extracting it cannot disagree with the printed
trace, and this package never carries a second copy of the answer (D-034).
"""
