# Smaller starting order for the live Higgs-jet transport

The notebook now starts physical transport with guard 20/order 16 at the same
requested 20-digit accuracy. All sixteen transports take **68.106 seconds** in
the frozen native solver experiment, compared with 91.512 seconds for guard
30/order 32 and 270.213 seconds for the original guard 60/order 96. The independent
uncertainty, truncation and precision checks are unchanged. Automatic boundary
generation keeps its original settings.

| Initial guard/order | All sixteen transports | Final bits/order | Master accuracy estimate | EW square / interference / HEFT square estimates |
| --- | ---: | ---: | ---: | --- |
| 60 / 96 | 270.213 s | 349 / 128 | 37 digits | 35 / 36 / 47 digits |
| 30 / 32 | 91.512 s | 249 / 64 | 37 digits | 35 / 36 / 47 digits |
| **20 / 16, selected** | **68.106 s** | **216 / 48** | **28 digits** | **35 / 36 / 47 digits** |
| 10 / 16, exploratory | 64.383 s | 183 / 48 | 20 digits | 30 / 31 / 47 digits |

Guard 10 passed the numerical checks too. Its further 5.5% time reduction uses
the remaining master-accuracy margin, so the notebook selects guard 20. These
are achieved estimates from the numerical checks, not claims based on working
precision. The independent EW-square reference still supports 19 comparison
digits.

Two separately written exact-rational checkers compare all 4,360 final complex
coefficients, all 8,720 inserted complex coefficients, eight form factors and
three observables with the original high-order run. Every difference fits the
combined admitted errors. Each run's individual errors and differences meet
the requested mixed 20-digit tolerance for master integrals and relative
20-digit tolerance for the nonzero factors and observables. All sixteen cases
take one accepted step with no rejected steps; the cache grows from 16 to 48
entries.

Original supplied values, 40-digit input caps, coordinates, root sheets, epsilon
ranges and provenance are unchanged. Derived cached boundaries inherit their
own achieved cap, which is 28 in the selected run versus 37 previously. The
independent checker permits that metadata difference only when it equals each
derived entry's achieved cap and remains at least 20; it does not weaken the
numerical comparisons. The original checker initially stopped at that
metadata-equality assertion; both updated checkers retain the numerical
inequalities and record the differing derived caps.

The native notebook controller also passed from an empty cache using the
shipped boundaries: all 4,360 reference comparisons, eight form factors, three
observables, exact binary reload and warm reuse. Boundary generation was
replaced by a failing hook, and comparison references were read only after the
calculation. Its measured stages were:

| Controller stage | Time |
| --- | ---: |
| Supplied-boundary import | 0.626 s |
| Initial transport, including checkpoint I/O | 75.295 s |
| First amplitude stage | 37.175 s |
| Binary reload and repeated transport | 2.354 s |
| Warm transport | 0.498 s |
| Warm amplitude | 2.668 s |

This controller uses the same release development extension as the earlier
supplied-notebook acceptance, pinned to RustFlow `7096ba8`. Its amplitude
constructor predates the separate staged-contraction optimization. The solver
experiment uses the frozen executable from the earlier compact-profile report.
Neither is a benchmark of the newer portable build or a new distribution wheel.
The updated controller, bundle and exporter suite passes 75 tests; Marimo lint
and whitespace checks pass.

A fresh process then copied the completed bank and evaluated the nearby point
`s + 1/100000`. All sixteen configurations selected the previous physical
destinations, retained 25 verified digits from 28-digit cached sources, and
inserted 32 points (48 to 80 cache entries). Nearby transport took 64.189 seconds
including checkpoint I/O. The newly projected observables retained 23/24/47
relative-digit estimates. This check tests source reuse and propagated accuracy;
it does not compare the changed physical point with the original-point
reference. Its first amplitude stage took 37.024 seconds with the same older
constructor. The archived original bank remained unchanged.

These are single observations per profile on one CPU of the shared EPYC 9754
host. The original baseline used CPU 43; both new profiles and the notebook
controller used CPU 44. Solver times exclude setup, Python, checkpoint I/O,
browser downloads and rendering. Controller times include checkpoint I/O but
exclude interpreter/model setup. They are **native measurements, not WASM
timings**. No change to Symbolica, HEPKit, RustRed, arithmetic algorithms or
thread counts was required for this profile change.

`report.json` records the artifact identities, settings, timings and checks.
`evidence.tar.gz` contains the frozen numerical outputs, exact-rational checkers,
controller sources and reports, and hashes of every archived file. It excludes
the executables and licensed configuration. The source-sensitive binary caches
retain their existing runtime compatibility checks.
