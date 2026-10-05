# Exact division rejected as statistically inconclusive

The physical nonplanar Higgs–jet boundary calculation exposed a Symbolica
`AtomField` bug at commit `942bd2c0cd2ef16414d69c176fc9eff7b21c0ab2`.
`try_div` and `try_inv` call `SelfRing::is_zero` on the denominator, which uses a
statistical test even when the field sets `statistical_zero_test: false`.
An inconclusive statistical result is consequently treated as zero. The
failure does not require RustFlow, floating-point inputs, or a singular matrix.

The [standalone example](../examples/symbolica_exact_division.rs) uses only
Symbolica. Its denominator is exactly `x + 1 - y`, represented as
`x + (10^40 - 10^40*y)/10^40`. The configured exact zero check correctly returns
false, but the unpatched division and inversion return `None`. A diagonal
matrix with entries `[denominator, 1, 1, 1]` then panics inside the Bareiss
determinant with “Inexact division”. Four is the smallest matrix dimension
using that algorithm; the message here denotes a rejected denominator, not
an actual nonzero remainder in polynomial division.

```sh
cargo run --release --locked --example symbolica_exact_division -- scalar
cargo run --release --locked --example symbolica_exact_division -- matrix
cargo run --release --locked --example symbolica_exact_division -- normalized
cargo run --release --locked --example symbolica_exact_division -- polynomial
```

To build independently of RustFlow and RustRed, copy that one source file into
a temporary Cargo project with only the pinned Symbolica dependency:

```sh
mre_dir=$(mktemp -d)
cp examples/symbolica_exact_division.rs "$mre_dir/repro.rs"
cat > "$mre_dir/Cargo.toml" <<'TOML'
[package]
name = "symbolica-exact-division-mre"
version = "0.0.0"
edition = "2024"

[[bin]]
name = "repro"
path = "repro.rs"

[dependencies]
symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "942bd2c0cd2ef16414d69c176fc9eff7b21c0ab2", default-features = false, features = ["integer-gmp", "float-mpfr"] }

[patch.crates-io]
numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "942bd2c0cd2ef16414d69c176fc9eff7b21c0ab2" }
graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "942bd2c0cd2ef16414d69c176fc9eff7b21c0ab2" }
TOML
cargo run --manifest-path "$mre_dir/Cargo.toml" -- scalar
cargo run --manifest-path "$mre_dir/Cargo.toml" -- matrix
```

On the unpatched owner, `scalar` and `matrix` fail. The `normalized` and
`polynomial` controls pass: explicitly expanding the input or using Symbolica's existing
`RationalPolynomialField` matrix gives the exact determinant `1 + x - y`.
Expansion is a workaround for this example, not a documented requirement for
exact `AtomField` division. `Matrix::det_in_place` also avoids the failure for
the original indicial matrix, but its elimination path does not retain row-swap
parity, so it is not a general determinant replacement.

The [owner patch](../scripts/patches/symbolica-exact-division.patch) changes both
admission checks to `<Self as Ring>::is_zero(self, ...)`, respecting the field's
chosen zero-test policy. It includes regressions for division/inversion,
rejection of literal zero, and the four-dimensional diagonal determinant.
Apply it to an isolated Symbolica checkout, and keep Symbolica, Numerica and
Graphica resolved from that same checkout when validating a consumer.
The isolated owner fix is committed locally as
`1fbdb0a92dc40e791a810e55552951a7b57166ad`, based on `942bd2c`; both selected
`domains::atom::exact_division_tests` passed. Its unit-test build enables
`native_code_generation`, as required by existing unrelated Symbolica test
imports. No cached upstream checkout was modified.

| Isolated check | Unpatched `942bd2c` | Patched owner |
|---|---|---|
| Exact scalar division and inverse | Both return `None` | Both pass exact multiplication checks |
| Unexpanded diagonal matrix | Bareiss panic | Exact `1 + x - y` |
| Expanded matrix control | Pass | Pass |
| Rational-polynomial matrix control | Pass | Pass |
| Original nonplanar indicial matrix | Bareiss panic | Expected exact factorization |

The original nonplanar four-dimensional indicial matrix was also reproduced
independently with Symbolica alone. With the owner patch it returns exactly
`(lambda - 2*eps)^3 * (lambda - 2*eps + 1/2)`. No changes to RustFlow's matrix
algorithm or expression normalization are needed for this correction.

Read-only upstream verification on 2026-10-05 found main at `75f8350094b90254ee71dc2a391fde0d14b0204a`
and community at `98794d0d7337ba2b08e4c046dde584ad7fc1ce10`. Both contain the
same division checks and matrix determinant implementation. The two later
community commits only affect Python NumPy handling and the license-server
address.

Local diagnostic artifacts are in `target/mre/symbolica-bareiss/`:
`minimal-old-*.log` records the failing cases and passing controls;
`minimal-patched-*.log` records the passing patched modes,
`patched-owner-tests-native.log` records the two owner tests, and
`patched-original-matrix.log` records the exact nonplanar determinant. The
original failure is retained in `target/gg-hg-nonplanar-physical-native.log`.
These small checks establish the owner fix; they do not certify the complete
nonplanar boundary fit or the final amplitude.
