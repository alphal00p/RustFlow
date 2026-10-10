The exact unimodular rerouting check passes, but the selected occupied-cut preparation fails at the frontier limit. It produces no numerical prediction.

The control represents original cut `[3]` as new cut `[0]`, retaining both original targets. The loop map has determinant −1 and absolute Jacobian one; incidence, chemical assignments, denominator powers and the fixed original raised numerator are preserved. Original raised edge 0 becomes edge 3. No extra cut or Wick phase arises. See `routing-audit.json` for the exact polynomial checks.

Using the unchanged captured K executable, the preparation reaches frontiers `20, 124, 603, 984, 1025` over five rounds, leaving 1,054 unresolved terms. Harness time was 65.61 s; measured process wall time was 65.63 s with peak child RSS 90,832 KiB. Compilation and Nix setup are excluded. This was a selected-sector, single-profile diagnostic with a 600 s limit. Search settings match the optional K schedule; the input representation, selected test/cut, output path and profile selection differ from K's full-amplitude request.

Routing, physical-edge order and numerator completion coordinates change together. The failure establishes neither numerical equivalence nor a benefit from any isolated representation choice. There is no full three-loop acceptance, imported rule corpus or performance comparison.

`summary.json` binds the outcome, exact input and captured build. `archive-map.json` records deterministic gzip archives with original/compressed hashes and byte counts. Inputs, settings, bindings, failure, and round 000/004 provisional proof pairs and their linked discovery records remain raw for subsequent source-only diagnostics. All archives are verified by exact restoration. `artifact-manifest.json` freezes this report.
