{
  description = "symbolica-amflow development environment";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  outputs = { self, nixpkgs }: let
    systems = [ "x86_64-linux" "aarch64-linux" ];
  in {
    devShells = nixpkgs.lib.genAttrs systems (system: let
      pkgs = import nixpkgs { inherit system; };
    in { default = pkgs.mkShell {
      packages = with pkgs; [ rustc cargo rustfmt clippy gcc gnum4 gnumake pkg-config perl git python3 ripgrep ];
      SYMBOLICA_HIDE_BANNER = "1";
    }; });
  };
}
