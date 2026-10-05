# Building the community extension from published inputs

This setup consumes RustFlow from the community manifest's published Git pin
and reconstructs the required HEPKit changes from public sources. Activate the
community repository's native Rust/Python build toolchain; the setup snippets
use Python 3.11 or newer. Run them from a fresh community integration checkout,
using a new directory for the build inputs. No particular sibling layout is
required.

The runtime pin is `4c48358185ec68de86f6acc6b195a201adbe2bae`. The extra RustFlow
checkout below supplies patch files and the Cargo configuration template; Cargo
loads the runtime itself from the community host's Git dependency. The locked
RustRed revision is `7c1ed03722b8c05daf60c89ba4ecc79457ed2ada`.

```sh
integration_deps="$PWD/../loop-integration-dependencies"
mkdir -p "$integration_deps"
integration_deps=$(cd "$integration_deps" && pwd)

git clone https://github.com/alphal00p/RustFlow "$integration_deps/RustFlow-inputs"
git -C "$integration_deps/RustFlow-inputs" checkout --detach 4c48358185ec68de86f6acc6b195a201adbe2bae
git clone https://github.com/alphal00p/gammaloop "$integration_deps/HEPKit"
git -C "$integration_deps/HEPKit" checkout --detach 9d086cec7971005ec7244b43e8fcea2a403c18ef

for patch in \
    hepkit-literal-substitution.patch \
    hepkit-isotropic-rank32.patch \
    hepkit-parameter-registration.patch \
    hepkit-linnet-strongly-connected.patch \
    hepkit-python-borrow-access.patch \
    hepkit-exact-parameters.patch
do
    patch_file="$integration_deps/RustFlow-inputs/scripts/patches/$patch"
    git -C "$integration_deps/HEPKit" apply --check "$patch_file" || exit
    git -C "$integration_deps/HEPKit" apply "$patch_file" || exit
done
```

Use a separate Cargo configuration for this build. The full patch table keeps
FeynKit, Linnet, Spenso, Idenso and their companion crates on one native owner.
The script generates that table from the published template, quoting the chosen
paths. It creates no RustFlow or Symbolica source override.

```sh
export CARGO_HOME="$integration_deps/cargo-home"
mkdir -p "$CARGO_HOME"
python3 - "$integration_deps/RustFlow-inputs" "$integration_deps/HEPKit" "$CARGO_HOME/config.toml" <<'PY'
import json
import pathlib
import sys
import tomllib

inputs, owner, destination = map(pathlib.Path, sys.argv[1:])
template = tomllib.loads((inputs / ".cargo/config.example.toml").read_text())
source = "https://github.com/alphal00p/gammaloop"
lines = [f"[patch.{json.dumps(source)}]"]
for name, entry in template["patch"][source].items():
    path = entry["path"].replace("/path/to/patched-hepkit", str(owner.resolve()))
    lines.append(f"{name} = {{ path = {json.dumps(path)} }}")
with destination.open("x") as output:
    output.write("\n".join(lines) + "\n")
PY
```

The checked-in community lock records the development path entry for RustFlow.
Resolve once without `--locked` to replace it with the published Git source;
subsequent operations use the resulting lock. This requires neither a local
RustFlow runtime override nor access to an unpublished owner commit.

```sh
cargo metadata --format-version 1 > "$integration_deps/resolved-graph.json"
cargo metadata --locked --offline --format-version 1 > "$integration_deps/locked-graph.json"
python3 - "$integration_deps/locked-graph.json" <<'PY'
import json
import sys

metadata = json.load(open(sys.argv[1]))
active = {node["id"] for node in metadata["resolve"]["nodes"]}
for name, revision in [
    ("symbolica-amflow", "4c48358185ec68de86f6acc6b195a201adbe2bae"),
    ("rustred-core", "7c1ed03722b8c05daf60c89ba4ecc79457ed2ada"),
    ("symbolica", "c3408e4ba1d3bdd4ea55678fad50e27009be13d4"),
    ("numerica", "c3408e4ba1d3bdd4ea55678fad50e27009be13d4"),
    ("graphica", "c3408e4ba1d3bdd4ea55678fad50e27009be13d4"),
]:
    matches = [p for p in metadata["packages"] if p["name"] == name and p["id"] in active]
    assert len(matches) == 1, (name, "expected one owner")
    source = matches[0]["source"]
    assert source and source.endswith("#" + revision), (name, source)
    print(name, source)
PY
```

With the community Python environment activated, build with
`maturin develop --release --locked`, then run the installed integration tests.
The host's tracked `.cargo/config.toml` already declares
`RUSTFLOW_WORKSPACE_MANIFEST`; preserve it so the shared extension fingerprints
the actual resolved graph. Keep the generated Cargo configuration and lock
with this build environment for later rebuilds and cache restarts.

On 2026-10-05, all six patches above were applied from published RustFlow
`6e86afd` to a disposable archive of public HEPKit `9d086ce`. All twelve affected
files matched the validated local owner `fc9ee6aa5` byte for byte. This verifies
the public patch reconstruction; a full clean-machine native build remains a
separate validation step. The local audit is recorded in
`target/published-hepkit-patch-audit.json`.

The HEPKit patch files remain required until those owner changes are available
upstream. Symbolica uses its official fixed community revision. Source-sensitive
cache keys include the resolved graph and configuration, so caches from another
development graph retain separate provenance.
