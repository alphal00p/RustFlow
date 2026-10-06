# Portable Higgs-plus-jet starting values

`HiggsJetIntegralSystem.kinematic_transport(options=None)` returns the existing generic `KinematicTransport` interface backed by the same canonical connection. Its `add_boundary` method accepts native arbitrary-precision coefficients and absolute errors, exact coordinates, a root germ, an accuracy cap and provenance. Importing supplied values performs the ordinary point, dimension, range, root and evidence validation. The fixed system requires `D = 4 - 2 epsilon` and `+i0`.

`HiggsJetIntegralSystem.mathematical_fingerprint` identifies the immutable scientific inputs independently of dependency versions, runtime builds, numerical backends and symbol namespaces. A portable bundle must compare this value before importing its values. The certificate identifies the mathematics; the numerical accuracy remains the supplied source's evidence. Normal binary caches retain their source-sensitive identities and compatibility checks.

The fingerprint is lowercase hexadecimal BLAKE3. Every field below is encoded as its byte length, an unsigned 64-bit little-endian integer, followed by its bytes. Text fields are UTF-8. File fields use the exact embedded bytes, without JSON parsing or newline conversion:

1. `symbolica-hep-integration:higgs-jet-mathematics:v1`
2. `Planar_EW1` or `NP_EW1`
3. `D=4-2*eps;+i0;measure=exp(2*eps*EulerGamma);mV2=mu2=1;basis=ordered-plugin-canonical;root-sheets=principal-or-opposite`
4. `integral-systems.json`, followed by that file from `fixtures/gg-hg`
5. `plugin-physical-map.json`, followed by that file from `fixtures/gg-hg`
6. `physical-configurations.json`, followed by that file from `fixtures/gg-hg`

Each filename and its contents are separate length-prefixed fields. A mathematical change to normalization, basis interpretation or branch conventions must update the convention tag or schema even if the data files do not change. The configuration file makes the certificate conservative when starting points or their declared sheets change.

Portable numerical encodings should retain exact integer ratios and recorded binary precision. Import them through Symbolica's native `Float` and `ComplexFloat` objects with their original error estimates and accuracy caps. Never convert them through Python `float` or `complex`.

The adapter introduces no transport, reduction, boundary generation, thread creation or new numeric owner. It prepares a shared handle for supplied-value import and subsequent numerical continuation. It does not itself establish a browser build; that remains the host's responsibility.
