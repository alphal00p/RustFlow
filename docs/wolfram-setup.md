# Wolfram runtime on the Bern cluster

Wolfram 15.0.1 is installed and successfully obtains a network license from
`itplic.itp.unibe.ch`. No Bern ID, Wolfram account login, or email address was
needed. Both the kernel and `wolframscript` have been exercised, including
arbitrary-precision arithmetic. The license remains dependent on access to the
Bern server during use.

The user-local launchers are:

```sh
/home/ben/.local/bin/bern-wolframscript -code '2+2'
/home/ben/.local/bin/bern-wolfram -noinit -noprompt -script /absolute/path/test.wls
```

These commands work without entering the Rust development environment. They
select the installed runtime and export `WOLFRAMINIT` so that child kernels,
including AMFlow's solver processes, inherit the same license configuration.
The configuration is `/home/ben/.config/symbolica-amflow/mathpass`, containing
only `!itplic.itp.unibe.ch`. It is readable only by its owner. The persistent Nix
GC root is `/home/ben/.local/share/symbolica-amflow/wolfram-15.0.1`.

The `wolframscript` launcher selects the patched native kernel executable
explicitly. Passing Nix's shell wrapper as `-local` caused WolframScript to try
to parse that wrapper as Wolfram Language; the native executable avoids this.

## Installation provenance and reproduction

The installer came from the official [Wolfram Download Center](https://www.wolfram.com/download-center/),
using the Linux download without local documentation. Its filename is
`Wolfram_15.0.1_LIN.sh`, size 2,644,745,046 bytes, and SHA-256 is
`ecde452688f481318dc18dd0fbc8491998becdee88451327c2eea72c4bf7705e`.
The proprietary installer and runtime are not committed to this repository.

The standard Linux installer assumes paths such as `/bin/bash` that this NixOS
host does not provide. The repository's pinned Nixpkgs package installs the
runtime and adjusts executable/library paths without changing its licensing.
CUDA support is disabled. To rebuild from a downloaded installer:

```sh
nix build --impure --file scripts/wolfram-runtime.nix \
  --argstr installer /absolute/path/Wolfram_15.0.1_LIN.sh \
  --out-link target/wolfram-runtime
```

The [runtime validation report](../reports/validation/2026-10-04-wolfram-runtime.json)
records the exact installer, Nixpkgs revision, installed store hash, license
connection, and successful kernel/CLI checks. Wolfram documents the network
server entry in its [mathpass guide](https://support.wolfram.com/112) and the
`-pwfile`/`WOLFRAMINIT` settings in its [kernel reference](https://reference.wolfram.com/language/ref/program/WolframKernel.html.en).

## AMFlow as an oracle

The licensed runtime has run the pinned original AMFlow 2.0 package with Kira
3.1 and Fermat 7.9b. Automatic Euclidean bubble and two-loop sunset calculations
passed 20-digit analytic comparisons. See [oracle usage and results](upstream-oracle.md)
for exact inputs, commands, live Rust comparisons, and remaining coverage limits.
These tests exercise the Mathematica implementation and its automatic reduction
and boundary workflow; the earlier standalone C++ comparisons remain separate.
