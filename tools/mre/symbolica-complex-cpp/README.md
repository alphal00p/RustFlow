# Scalar complex constants generate invalid C++

This reproducer uses only the public scalar Symbolica Python API. It fails on
official Symbolica community commit
[`c3408e4ba1d3bdd4ea55678fad50e27009be13d4`](https://github.com/symbolica-dev/symbolica/commit/c3408e4ba1d3bdd4ea55678fad50e27009be13d4).
No HEPKit, Spenso, RustFlow, numerical reference data, or local Symbolica patch is
needed. A C++ compiler and an ordinary licensed Symbolica installation are
required; no license key is included here.

In a Python environment built from that revision, run:

```sh
python tools/mre/symbolica-complex-cpp/reproduce.py
```

To build an independent Python installation from the pinned source, the usual
Rust and Symbolica native build prerequisites must be installed first:

```sh
python -m venv /tmp/symbolica-complex-cpp-venv
. /tmp/symbolica-complex-cpp-venv/bin/activate
python -m pip install numpy 'git+https://github.com/symbolica-dev/symbolica.git@c3408e4ba1d3bdd4ea55678fad50e27009be13d4'
python tools/mre/symbolica-complex-cpp/reproduce.py
```

The script evaluates `i*x` with `params=[x]`, then compiles the same evaluator
with the documented `number_type="complex"`, `inline_asm="none"`, and
`native=False` options. It exercises both `direct_translation=True` and
`False`. The imaginary unit is a constant, not an input parameter. For inputs
`2+3i`, `-4`, and `-5i`, both the interpreter and compiled library must return
`-3+2i`, `-4i`, and `5` exactly. The interpreter control passes, but both C++
compilations fail on the recorded revision; the script exits with status 1.
It exits 0 only if both compiled evaluations also satisfy these checks.

The generated code contains:

```cpp
template<typename T>
void complex_direct_complexf64_gen(T* params, T* Z, T* out) {
    Z[2] = params[0]*T(T(0e0), T(1e0));
    out[0] = Z[2];
}
```

The wrapper instantiates `T=std::complex<double>`, whose two-component
constructor requires two real components, not two complex values. The compiler
therefore reports no matching constructor for
`std::complex<double>::complex(std::complex<double>, std::complex<double>)`.
Symbolica's `src/evaluate/export.rs` implementation of
`Complex<T>::export_wrapped_with` wraps each component with the outer complex
wrapper before constructing the complex constant. The public Python compile
method correctly dispatches to `export_cpp::<Complex<f64>>`.

By default all generated files use a temporary directory. To retain them:

```sh
python tools/mre/symbolica-complex-cpp/reproduce.py \
  --output-directory /tmp/symbolica-complex-cpp-output
```

`provenance.json` records the tested revision, installed extension hash, and
observed failure. This is a C++ export issue. RustFlow's native arbitrary-
precision MPFR loop solver does not use this exporter and is unaffected by
this particular failure. The retained upstream tensor compilation test remains
a failing compatibility gate until the exporter is repaired.
