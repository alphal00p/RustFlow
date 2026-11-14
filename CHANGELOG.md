# Changelog

## Unreleased

### Added
- Dispersive mode (`src/dispersive`, `examples/dispersive.rs`): values at fixed `D` from the discontinuity of one line's
  auxiliary-mass equation (exact e-block from a RustRed DE; mod-p structure and MPFR Frobenius families at zero and infinity; soft
  boundary from lower-loop samples with a full-rank and consistency gate; optional native four-loop prism-lane stage; matching at `S`
  and dispersion integral). Every stage checkpoints and resumes (any thread count; exit code 75 on SIGTERM); two-precision per-entry
  digits; per-stage cost records; `--progress` from the checkpoint files. See `docs/dispersive.md`.
- Fixed-precision patch kernel (`src/patch_kernel`, `examples/patch_transport.rs`): planned Taylor-patch transport of exact rational
  systems with exact shifts, bit-identical row-parallel recurrence, exact checkpoints with fingerprints, overwrite guard and lock,
  fail-fast on non-finite values or divergent patches; `TransportConfig::dispersive` and `::generic`; `RationalSystem` from entries,
  atoms, `DifferentialSystem` or `CompiledSystem::patch_kernel_system`.
- Validation at a five-loop vacuum node (class 16359, eps = 1/150000, 569 masters): 59.2-63.8 digits agreement with the reference
  tool's sample on the five target masters, 0.57x its CPU time; thread-count and kill/resume identity controls.
