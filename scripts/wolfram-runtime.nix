# Build the official Linux installer with the repository's pinned Nixpkgs.
# Supply an installer downloaded from Wolfram; no license is embedded here.
{ installer, version ? "15.0.1" }:
let
  lock = builtins.fromJSON (builtins.readFile ../flake.lock);
  nixpkgs = builtins.getFlake "github:NixOS/nixpkgs/${lock.nodes.nixpkgs.locked.rev}";
  pkgs = import nixpkgs {
    config.allowUnfree = true;
  };
in
pkgs.mathematica.override {
  cudaSupport = false;
  versionInfo = { inherit version; lang = "en"; };
  source = builtins.path {
    path = installer;
    name = builtins.baseNameOf installer;
  };
}
