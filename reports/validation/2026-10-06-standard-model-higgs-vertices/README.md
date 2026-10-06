# Native Standard Model form factors and numerical-method citations

The effective ggHg declarations now extend HEPKit's native Standard Model through
its model interchange. Original SM data is preserved, declaration collisions
fail, and particle selection uses PDG identities. The declarations contain
symbolic tensors/couplings and zero placeholders, not numerical loop results.
Native scientific/methodological citations are cumulative and activate only
after successful preparation; the community host can merge the returned native
Citation records into Symbolica's existing registry.

Validation passed: two native amplitude tests, four Python/stub tests, focused
strict release Clippy, formatting/diff checks and exact asset regeneration.
The three Standard Model scalar kernels equal the legacy kernels exactly.
Python tests cover native model ownership, rejected collisions, six isolated
citation scenarios and exact supplied-boundary restart. One independent
community-bundle acceptance test remains explicitly ignored in this gate.

The generated integration stub uses the host's compatibility normalizer and
passes Python3.9 syntax parsing. Its generator is archived; the temporary Cargo
example was removed. Production changes and fixture hashes remain those frozen
before the native/Python gates.

See [report.json](report.json) and [evidence.tar.gz](evidence.tar.gz) for commands,
source/asset/stub hashes and accepted logs. Earlier development failures and the
old-runtime asset-only validation are preserved separately. New shared-host
native/browser notebook acceptance is a subsequent host gate; this milestone
does not claim the existing browser wheel already includes the new APIs.
