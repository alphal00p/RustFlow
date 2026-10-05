# Optional native amplitude validation

This independent Cargo workspace keeps FeynKit generation and amplitude dependencies
outside the RustFlow library. It reproduces the native `gg → Hg` tensor contraction
and evaluates the recorded coherent W/Z form factors. Native evaluation needs no
Python; the optional `ufo-import` feature uses HEPKit's existing `UfoLoader` for
input conversion. The dependency layout expects `../hepkit` beside this checkout.

The external UFO must come from `mg5_higgs_ew_plugin` at commit
`89f64b93d0bdbd8ee90eb85737021189229b034a`. No upstream model or generated MG5
computational source is bundled here. The preparation helper checks the Git
revision and every tracked UFO file against that revision before copying. It
rejects modified/additional input files, symlinks, and existing output directories;
Python bytecode caches are ignored. Original model files are never edited or
imported by the helper.

## Prepare and scope the external input

Run from the RustFlow checkout, using a fresh work directory:

```sh
UFO_SOURCE=/absolute/path/to/mg5_higgs_ew_plugin/UFO_model_gggH
UFO_WORK="$PWD/target/hepkit-amplitude-input"
python3 tools/hepkit-amplitude/scripts/prepare_ufo.py prepare \
  --source "$UFO_SOURCE" --work "$UFO_WORK"
cargo test --manifest-path tools/hepkit-amplitude/Cargo.toml --release --bin scope-lorentz
cargo run --manifest-path tools/hepkit-amplitude/Cargo.toml --release \
  --bin scope-lorentz -- "$UFO_WORK"
python3 tools/hepkit-amplitude/scripts/prepare_ufo.py apply-scopes --work "$UFO_WORK"
```

`prepare` makes `UFO_model_gggH_python3`, `original-lorentz.json` and
`preparation.json`. Its recorded compatibility edits are Python 3 dictionary/print
syntax, an explicit function-library import, function-call whitespace, the native
ghost spin signature, and the four-scalar `V_1001` Lorentz slot metadata. The two
scalar structures are checked to be exactly `1`. All 45 tensor structure strings
must remain byte-for-byte identical at this stage. The script parses Python AST
literals; it does not execute the model or interpret tensor expressions.

The native `scope-lorentz` binary parses the original tensor strings with
Symbolica's token parser and preserves reciprocal momentum scopes with the native
`spenso::bracket` representation. It writes `scoped-lorentz.json` and
`scope-counts.json`. `apply-scopes` then creates `UFO_model_gggH_scoped`, replacing
only the named `Lorentz(..., structure=...)` string literals by their UTF-8 AST
spans. It rejects unknown/duplicate names or a changed prepared input and records
before/after SHA256 hashes for every file in `scoping-provenance.json`. Only
`lorentz.py` changes; original and prepared copies remain available for inspection.
Neither helper implements tensor contraction or changes the original tensor algebra.

## Import with the native loader

Use one consistent Python 3.13 environment for the optional importer. The validated
Python packages were Symbolica 3.0.0 and the pinned UFO loader below:

```sh
python3.13 -m venv "$PWD/target/hepkit-amplitude-python"
UFO_PYTHON="$PWD/target/hepkit-amplitude-python/bin/python"
env -u PYTHONPATH "$UFO_PYTHON" -m pip install \
  'symbolica==3.0.0' \
  'git+https://github.com/alphal00p/ufo_model_loader.git@70ddee6b416f8c8b340e0d087646d77095c5d24b'
unset PYTHONHOME PYTHONPATH PYO3_CONFIG_FILE
export PYO3_PYTHON="$UFO_PYTHON"
export PYTHONPATH="$("$UFO_PYTHON" -c 'import site; print(site.getsitepackages()[0])')"
cargo run --manifest-path tools/hepkit-amplitude/Cargo.toml --release \
  --features ufo-import --bin native-ufo-export -- \
  "$UFO_WORK/UFO_model_gggH_scoped" "$UFO_WORK/native-model.json"
```

In the validated Nix environment, `.dev-env` supplied Python 3.14 site-packages.
Combining those paths with the Python 3.13 embedded library crashed the prototype.
Replace that `PYTHONPATH` with the chosen environment's site-packages as shown;
the embedded interpreter does not automatically activate a virtual environment.
If a custom `PYO3_CONFIG_FILE`
is used, its executable/library must agree with `PYO3_PYTHON` and the runtime
Python 3.13 paths. This constraint concerns the optional Python-backed input
conversion; `native-amplitude` and `evaluate-me` operate on native model JSON.

`native-ufo-export` calls the existing HEPKit loader with restriction `full` and
model simplification disabled. It preserves the full UFO model and does not
select binary64 default parameters for the scientific evaluation.

## Reproduce the recorded scalar comparison

The durable report contains the exact phase point, normalized form-factor strings,
their absolute allowances, and the independent MG5/ALOHA binary128 output. Extract
those recorded inputs without numerical conversion:

```sh
python3 tools/hepkit-amplitude/scripts/prepare_ufo.py extract-validation-inputs \
  --report reports/validation/2026-10-04-gg-hg-squared-matrix-element.json \
  --output "$UFO_WORK/validation"
cargo run --manifest-path tools/hepkit-amplitude/Cargo.toml --release \
  --bin native-amplitude -- "$UFO_WORK/native-model.json" "$UFO_WORK/tensor-output"
cargo run --manifest-path tools/hepkit-amplitude/Cargo.toml --release \
  --bin evaluate-me -- "$UFO_WORK" \
  "$UFO_WORK/validation/amplitude-input.json" "$UFO_WORK/validation/quad-result.json"
```

The extractor maps the report's `normalized_form_factors` entries to the native
harness input keys, retaining decimal strings and uncertainties exactly. It also
exports the original oracle result and records the source report/output hashes.
It does not regenerate the independent oracle. To evaluate new scientific inputs,
provide the same explicit input schema and independently generated oracle data.

The observable is pure `|A_W+A_Z|²`, averaged over the two incoming gluons' spins
and colors and summed over the final state. It is not the HEFT–EW interference.
The native ownership chain is FeynKit process generation, amplitude squaring,
spin/color sums, Idenso contraction and SU(3) simplification, native kinematics,
model coupling expansion, and Symbolica MPFR evaluation. No replacement tensor
algebra or generated MG5 computational routines are included in this workspace.

The recorded source-error estimate supports 19 relative digits for this matrix
element. Its approximately `3.03e-33` relative agreement with the independent
binary128 contraction checks the implementation; it does not raise that physical
accuracy claim. Form-factor uncertainty and arithmetic refinement remain separate.


## HEFT LO and HEFT–EW interference

The optional `evaluate-heft` binary reuses `scalar.rs` with `evaluate-me`: native HEPKit generates the amplitude, performs state sums, and delegates tensor contraction to Idenso. Symbolica separates the resulting scalar squared amplitude by an auxiliary HEFT coefficient. Degrees zero, one, and two give the EW square, `2 Re(A_HEFT conjugate(A_W+A_Z))`, and the HEFT square. The exact scalar reconstruction is checked, including absence of hidden marker dependence. No additional factor of two or incoming-state average is applied.

This benchmark uses the plugin's **infinite-top HEFT** LO amplitude. It does not validate finite-top QCD. The exact HEFT scalar inputs are prepared from the externally supplied pinned C++ coefficient bridge; the original UFO's numerical placeholder values are never used. `prepare_heft.py` verifies the external source hash, reads its four scalar assignments with a closed arithmetic AST, and evaluates only exact rational operations at the supplied phase point. It neither imports upstream code nor implements tensor algebra.

After the existing preparation/import steps have produced `$WORK/native-model.json` and the coherent inputs, run from the repository root:

```sh
python tools/hepkit-amplitude/scripts/prepare_heft.py \
  --source "$PLUGIN/ComputationFormFacHEFT/fortran_bridge_HEFT_ggHg.cpp" \
  --coherent-input "$INPUTS/amplitude-input.json" \
  --oracle-report reports/validation/2026-10-05-gg-hg-heft-interference.json \
  --output "$HEFT_INPUTS"

cargo run --release --manifest-path tools/hepkit-amplitude/Cargo.toml \
  --bin evaluate-heft -- \
  "$WORK" "$INPUTS/amplitude-input.json" \
  "$HEFT_INPUTS/exact-heft.json" "$HEFT_INPUTS/heft-oracle.json" \
  "$HEFT_RESULT"
```

`$HEFT_INPUTS` must be a new directory. All phase coordinates, HEFT rational coefficients, oracle outputs, and EW inputs retain string representations until exact Symbolica/MPFR evaluation. The tool rejects incomplete/duplicate HEFT coefficients, a different phase point, finite-precision HEFT inputs, failed arithmetic refinement, or a failed original comparison. The argument paths are supplied explicitly; no in-repository private model or target-directory artifact is assumed by the tool.

The original selections used for independent reference generation were:

```text
generate g g > g h GGGHEFT^2==2 QCD^2==6 GGHEFT=0 GGHEW=0 GGGHEW=0 @1
add process g g > g h GGGHEFT^2==1 GGGHEW^2==1 QCD^2==6 GGHEFT=0 GGHEW=0 @2
```

The report preserves ordinary binary64 outputs and native MG5/ALOHA MP outputs (113-bit arithmetic). The MP adapter changes numerical kinds, native MP routine names, and common-block layouts while keeping generated color, helicity, amplitude, order-selection, and averaging expressions unchanged. The original scalar HEFT bridge was also rerun with a numerical-kind-only `__float128`/`powq` adaptation. The checked arithmetic comparison uses a conservative `1e-30` relative tolerance.

At the recorded coherent point the interference has a conditional relative allowance of approximately `8.3214e-21`, supporting 20 relative digits. This follows from the exact linear polynomial in 16 real EW form-factor components: the native gradient propagates their supplied complex allowances, and a `1e-65` arithmetic reserve covers the much smaller 256/384-bit change and imaginary cancellation residual. This is empirical input-error propagation, not an interval certificate. The HEFT LO check establishes arithmetic agreement and refinement; it does not derive a precision guarantee merely from working bits. The separate pure EW square still has its previously recorded 19-digit relative bound.
