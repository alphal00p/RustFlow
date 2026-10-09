# Native interface validation after checkpoint bf60718

The CLI now invokes native complete-amplitude preparation and evaluation through
`finite_density::assembly::PreparedDensityFlow`. A required sector failure is an
error. The shared serializer preserves decimal precision, ordered targets,
individual cut contributions, contour and boundary provenance, exact conditions,
strict support-zero proofs and Laurent accuracy metadata. Values use raw,
unscaled Euclidean loop measures.

The Python-feature all-tests type check passed in 15.42 seconds. The dependency
owner ran the three shared interface unit tests successfully. This is compiled
binding evidence; loading the function in an actual embedding Python host has
not been tested. The compiled source snapshot is `../tangent-source-hashes.json`.

| Native command or check | Actual result | Process wall seconds | Peak RSS KiB |
| --- | --- | --- | --- |
| Nine input regressions, including all 108 exact oracle basis mappings | 9/9 passed | 0.1484 | 12336 |
| Four normalization regressions | 4/4 passed | 0.3065 | 12336 |
| Schema-version-2 preparation | Passed; performed=false, separate numerical interfaces identified | 0.0065 | 12320 |
| Complete tadpole sample, mu=1 > m=1/2 | Passed | 0.0401 | 12360 |
| Complete tadpole sample, mu=1/4 < m=1/2 | Passed, with strict occupied support-zero proof | 0.1802 | 12320 |
| Zero-epsilon sample | Correctly rejected, exit 1 before preparation | 0.0045 | 12320 |
| Complete tadpole Laurent command | Passed; both targets report 12 verified digits | 0.7367 | 12344 |

Measurements use the prebuilt executable or test binary and exclude compilation
and Nix startup. Other validation and build processes were active; these are
recorded smoke-test resources, not isolated throughput benchmarks. Exact commands,
exit codes, CPU times and complete output are in the adjacent `*-resources.json`,
`*-regressions.json` and `.log` files.

The sample tests use D=3, epsilon=1/2, with scalar and raised scalar powers. Their
independent complete-amplitude references are

```text
I1 = -max(m,mu)/(4*pi),
I2 = 1/(8*pi*m) below threshold and 0 above threshold.
```

All four sample comparisons pass relative error 1e-12 for nonzero values and
absolute error 1e-25 for the zero value. Above threshold, the raised vacuum and
occupied terms cancel to a residual 3.7577758385489893e-46. Successful output
identifies `compact_polynomial_moments` as the occupied construction; below
threshold it retains the exact disjoint-support proof. A single sample correctly
reports `verified_digits: null`.

The Laurent CLI uses D=4-2 epsilon and requests orders -1 and 0. For a=m²,
r=sqrt(mu²-a), L=ln((mu+r)/sqrt(a)), the independent coefficients above threshold
are

```text
I1[-1] = -a/(16*pi²),
I1[0]  = a/(16*pi²)*(gamma_E-1-ln(4*pi/a))
         - (mu*r-a*L)/(8*pi²),
I2[-1] = 1/(16*pi²),
I2[0]  = (ln(4*pi/a)-gamma_E-2*L)/(16*pi²).
```

These follow from the Euclidean Gaussian vacuum Gamma function and the compact
radial integral. The raised expression is differentiated at fixed original
numerator and chemical potential. The comparison script obtains pi from Machin's
arctangent identity and gamma_E from an independent rational-Bernoulli
Euler–Maclaurin expansion using 100-digit Decimal arithmetic. It does not read the
oracle or call production normalization, boundary or integration functions.
All four Laurent coefficients pass the 1e-12 relative comparison; the largest
observed relative discrepancy is 3.9417847681096216e-42. Native output separately
retains its independent fit changes, working bits, sample counts and refinement
count. The comparison and all original decimal strings are saved in
`tadpole-independent-comparison.json`.

Reproduce the public commands from the repository root:

```sh
target/release/rustflow finite-density-sample examples/finite_density/massive_one_loop_tadpole.json 1/2
target/release/rustflow finite-density examples/finite_density/massive_one_loop_tadpole.json
python3 tools/finite_density/compare_tadpole_cli.py
```

This validates complete one-loop interface execution. It does not resolve the
above-threshold massive two-loop closure failure, establish the three mandatory
four-loop predictions, or count as numerical comparison of a supplied oracle
record. Those mandatory acceptance gates remain open.
