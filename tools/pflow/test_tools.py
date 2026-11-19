"""Tests of the authenticated pflow tools (GATE 1: input authentication) on synthetic inputs.
run: python3 -I -m unittest discover -s tools/pflow -p 'test_*.py'   (needs mpmath)"""
import json
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as Fr
from pathlib import Path

H = Path(__file__).resolve().parent
sys.path.insert(0, str(H))
from mpmath import mp, mpf, gamma, zeta  # noqa: E402
from pflowio import Refused, eval_exact, read_masters, read_sample, rat, seal, sha256_canonical, sha256_file, canonical  # noqa: E402
import fitgrid  # noqa: E402


def g11(eps):
    return gamma(1 + eps) * gamma(1 - eps) ** 2 / gamma(2 - 2 * eps)


def run_tool(name, *args):
    p = subprocess.run([sys.executable, "-I", "-B", str(H / name), *map(str, args)], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


class Exact(unittest.TestCase):
    def test_canonical_is_language_independent(self):
        """the same string as RustFlow's pflow::cli::canonical_json (Rust test canonical_serialisation_is_exact)"""
        v = {"b": 0.5, "a": [1, 2.0, 0.25, -3], "s": "x\"y", "n": None, "t": True}
        self.assertEqual(canonical(v), '{"a":[1,2,1/4,-3],"b":1/2,"n":null,"s":"x\\"y","t":true}')

    def test_typed_refusals_of_the_canonicaliser(self):
        """GATE 1 re-check 2 (c, d): wide integers, NaN / Infinity and lone surrogates are typed refusals, never a traceback"""
        for v, msg in [({"a": 1 << 70}, "beyond 64 bits"), ({"a": float("nan")}, "non-finite"), ({"a": float("inf")}, "non-finite"),
                       ({"a": "x\ud800y"}, "lone surrogate")]:
            with self.assertRaises(Refused) as c:
                canonical(v)
            self.assertIn(msg, str(c.exception))
        from pflowio import load_json
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "nan.json"
            p.write_text('{"a": NaN}')
            with self.assertRaises(Refused) as c:
                load_json(p)
            self.assertIn("non-finite", str(c.exception))

    def test_rational_arithmetic(self):
        d = Fr(4)
        self.assertEqual(eval_exact("(d+1)/3", {"d": d}), Fr(5, 3))
        self.assertEqual(eval_exact("d*(1/3)", {"d": d}), Fr(4, 3))  # no binary floating point
        self.assertEqual(eval_exact("(2*d^2-3)/(d-1)^2", {"d": Fr(7, 2)}), Fr(2 * 49 - 12, 25) * 1)
        with self.assertRaises(Refused):
            eval_exact("1.5*d", {"d": d})
        with self.assertRaises(Refused):
            rat(1.5)


def write(dirp, name, obj):
    p = Path(dirp) / name
    p.write_text(json.dumps(obj))
    return p


DEFS = {"m1": {"family": "P3L_BE", "index": [1, 1, 1, 1, 0, 1, 1, 1, 0], "relation": "M = eps * E_G(key)"}}
KEY = [[0, 2, "1"], [0, 3, "1"], [1, 2, "1"], [1, 3, "1"], [2, 3, "1"]]  # a two-loop key (L = 2)


def answer(d0_listed=None, coef="(d-3)/(d+1)"):
    x = {"id": "X1", "family": "T", "index": [1], "sign": -1,
         "terms": [{"basis": "G[(1,1),(-2,1)]", "coefficient": coef}, {"basis": "M:m1", "coefficient": "1/(d-2)"}]}
    if d0_listed:
        x["at_d0"] = {d0_listed[0]: {"terms": [{"coefficient_exact": d0_listed[1]}, {"coefficient_exact": "2/3"}]}}
    return {"schema": "num-support.req-answer.v1", "meaning": "test", "master_definitions": DEFS, "integrals": [x]}


class Converter(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.d = self.t.name
        self.keys = write(self.d, "keys.json", {"ids": {json.dumps(KEY): "X1"}})

    def tearDown(self):
        self.t.cleanup()

    def masters(self, d0, err="1e-50", sources=(), edit=None):
        from pflowio import digest
        m = seal({"schema": "pflow-masters.v1", "d0": d0, "master_definitions_sha256": digest(DEFS),
                  "values": {"m1": {"value": "2.5", "abs_err": err}}, "sources": list(sources)})
        if edit:
            edit(m)
        return write(self.d, f"m-{d0.replace('/', '_')}.json", m)

    def test_edited_master_file_and_stale_source(self):
        """review GATE 1 A1: an edited master value (seal kept) and a stale source digest are refused; an unreachable source is a
        recorded downgrade in the table"""
        a = write(self.d, "a.json", answer())
        out = Path(self.d) / "p.json"
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", out, a, "--masters", self.masters("7/2", edit=lambda m: m["values"]["m1"].update(value="2.6")))
        self.assertEqual(rc, 2)
        self.assertIn("edited", err)
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", out, a, "--masters", self.masters("7/2", sources=[{"file": str(a), "sha256": "2" * 64}]))
        self.assertEqual(rc, 2)
        self.assertIn("stale", err)
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", out, a, "--masters", self.masters("7/2", sources=[{"file": "/nonexistent/x.json", "sha256": "2" * 64}]))
        self.assertEqual(rc, 0, err)
        t = json.loads(out.read_text())
        self.assertIn("not reachable", " ".join(t["downgrades"]))
        from pflowio import content_sha256
        self.assertEqual(t["content_sha256"], content_sha256(t))

    def test_value_and_propagated_uncertainty(self):
        a = write(self.d, "a.json", answer())
        rc, out, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, "--masters", self.masters("7/2"))
        self.assertEqual(rc, 0, err)
        t = json.loads((Path(self.d) / "p.json").read_text())
        self.assertEqual((t["schema"], t["d0"], t["unverified"]), ("pflow-provided.v2", "7/2", []))
        mp.dps = 80
        d = mpf(7) / 2
        eps = (4 - d) / 2
        c1 = (d - 3) / (d + 1)
        c2 = 1 / (d - 2)
        want = -(c1 * gamma(1 + eps) * gamma(1 - 2 * eps) + c2 * mpf("2.5")) * g11(eps) ** 2
        v = t["values"][0]
        self.assertLess(abs(mpf(v["value"]) / want - 1), mpf(10) ** -70)
        e = mpf(v["abs_err"])
        prop = abs(c2) * mpf("1e-50") * abs(g11(eps)) ** 2
        self.assertTrue(prop <= e <= prop * 1.02, (e, prop))  # the master's error, propagated (rounded up to 3 digits)

    def test_master_file_of_another_dimension_is_refused(self):
        a = write(self.d, "a.json", answer())
        out = Path(self.d) / "p.json"
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", out, a, "--masters", self.masters("15/4"))
        self.assertEqual(rc, 2)
        self.assertIn("d0 = 15/4, required d0 = 7/2", err)
        self.assertFalse(out.exists(), "nothing may be written")

    def test_definitions_digest_and_answer_coefficients(self):
        a = write(self.d, "a.json", answer())
        m = json.loads(self.masters("7/2").read_text())
        m["master_definitions_sha256"] = "0" * 64
        mf = write(self.d, "mbad.json", seal(m))
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, "--masters", mf)
        self.assertEqual(rc, 2)
        self.assertIn("master_definitions_sha256", err)
        # the answer lists its own exact coefficient at this d0: a disagreement is refused
        bad = write(self.d, "b.json", answer(("7/2", "1/7")))
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", bad, "--masters", self.masters("7/2"))
        self.assertEqual(rc, 2)
        self.assertIn("the answer lists", err)
        ok = write(self.d, "c.json", answer(("7/2", "1/9")))
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", ok, "--masters", self.masters("7/2"))
        self.assertEqual(rc, 0, err)
        # two answer files with different master definitions
        other = answer()
        other["master_definitions"] = {"m1": dict(DEFS["m1"], relation="M = E_G(key)")}
        o = write(self.d, "o.json", other)
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, o, "--masters", self.masters("7/2"))
        self.assertEqual(rc, 2)
        self.assertIn("differ", err)

    def test_need_digits_and_removed_stamp(self):
        a = write(self.d, "a.json", answer())
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, "--masters", self.masters("7/2", "1e-20"), "--need-digits", "40")
        self.assertEqual(rc, 2)
        self.assertIn("--need-digits", err)
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, "--masters", self.masters("7/2"), "--master-digits", "90")
        self.assertEqual(rc, 2)

    def test_legacy_master_file(self):
        a = write(self.d, "a.json", answer())
        lg = write(self.d, "legacy.json", {"m1": "2.5", "d0": "7/2"})
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, "--masters", lg)
        self.assertEqual(rc, 2, "a legacy file without a readable source run is refused")
        src = write(self.d, "X.json", {"schema": "pflow-route-b.v1", "d0": "7/2", "config": {"w": 60}})
        lg = write(self.d, "legacy2.json", {"m1": "2.5", "d0": "7/2", "from": str(src)})
        rc, _, err = run_tool("answer2prov.py", self.keys, "7/2", Path(self.d) / "p.json", a, "--masters", lg)
        self.assertEqual(rc, 0, err)
        t = json.loads((Path(self.d) / "p.json").read_text())
        self.assertEqual(t["unverified"], ["X1"])
        src2 = write(self.d, "X2.json", {"schema": "pflow-route-b.v1", "d0": "15/4", "config": {"w": 60}})
        lg = write(self.d, "legacy3.json", {"m1": "2.5", "d0": "7/2", "from": str(src2)})
        with self.assertRaises(Refused):
            read_masters(lg, Fr(7, 2))


FRAME = {"normalisation": "E, P^2=1", "graph": {"edges": [[0, 1], [0, 1]], "ext": [0, 1], "loops": 1}, "e_slot": 0, "masters": [[1, 1]], "boundary": None}
FSHA = sha256_canonical(FRAME)


def sample(d0, value, err="1e-40", accepted=True, frame=FRAME, target=30.0, tp=5, mutate=None, **kw):
    """a sealed v2 sample (frame body and digest, two-precision pass); mutate() runs BEFORE the seal (a different producer output)"""
    s = {"schema": "pflow-route-b.v2", "accepted": accepted, "acceptance": {"accepted": accepted, "failures": [] if accepted else ["x"]},
         "accept_diagnostics": False, "d0": d0, "frame": json.loads(json.dumps(frame)), "frame_sha256": sha256_canonical(frame), "target_digits": target,
         "config": {"two_precision": tp}, "downgrades": [], "inputs": {},
         "values": [{"m": 1, "master": [1, 1], "value": value, "abs_err": err}], "extras": []}
    s.update(kw)
    if mutate:
        mutate(s)
    return seal(s)


class Samples(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.d = self.t.name

    def tearDown(self):
        self.t.cleanup()

    def test_refusals(self):
        ok = write(self.d, "ok.json", sample("7/2", "1.5"))
        r = read_sample(ok, d0="7/2", frame_sha=FSHA, need_target=20)
        self.assertEqual((r["values"][(1, (1, 1))][0], r["downgrades"]), ("1.5", []))
        other = dict(FRAME, e_slot=3)

        def edited(f):
            s = sample("7/2", "1.5")
            f(s)
            return s
        for name, s, kw, msg in [
            ("rej", sample("7/2", "1.5", accepted=False), {}, "REJECTED"),
            ("nov", sample("7/2", "1.5", mutate=lambda s: s.pop("accepted")), {}, "no acceptance verdict"),
            ("diag", sample("7/2", "1.5", accept_diagnostics=True), {}, "accept-diagnostics"),
            ("d0", sample("15/4", "1.5"), {"d0": "7/2"}, "required d0"),
            ("frame", sample("7/2", "1.5", frame=other), {"frame_sha": FSHA}, "frame_sha256"),
            ("tgt", sample("7/2", "1.5", target=12.0), {"need_target": 18}, "needed"),
            ("v1", {"schema": "pflow-route-b.v1", "accepted": True, "d0": "7/2", "values": [], "config": {"w": 40}}, {}, "schema"),
            # review GATE 1 A1: content edits after the producer's seal
            ("edited value", edited(lambda s: s["values"][0].update(value="1.6")), {}, "edited"),
            ("edited abs_err", edited(lambda s: s["values"][0].update(abs_err="1e-60")), {}, "edited"),
            ("edited frame", edited(lambda s: s["frame"].update(e_slot=3)), {}, "edited"),
            ("unsealed", edited(lambda s: s.pop("content_sha256")), {}, "unsealed"),
            # a frame body that does not match its digest even under a fresh seal
            ("frame body", sample("7/2", "1.5", mutate=lambda s: s["frame"].update(e_slot=3)), {}, "frame was edited"),
            ("stale input", sample("7/2", "1.5", inputs={"de": [{"file": __file__, "sha256": "3" * 64}]}), {}, "stale"),
        ]:
            with self.assertRaises(Refused) as c:
                read_sample(write(self.d, name.replace(" ", "_") + ".json", s), **kw)
            self.assertIn(msg, str(c.exception), name)
        # downgrades are reported, not hidden: no two-precision pass, an unreachable input
        r = read_sample(write(self.d, "down.json", sample("7/2", "1.5", tp=0, inputs={"de": [{"file": "/nonexistent", "sha256": "4" * 64}]})))
        self.assertEqual(len(r["downgrades"]), 2, r["downgrades"])
        # the legacy path admits ONLY pinned archived four-loop nodes: any other v1 file is refused
        v1 = write(self.d, "v1-new.json", {"schema": "pflow-route-b.v1", "accepted": True, "d0": "7/2", "values": [], "config": {"w": 40}})
        with self.assertRaises(Refused) as c:
            read_sample(v1, legacy_pins={sha256_file(v1)})
        self.assertIn("only for the 12 archived four-loop nodes", str(c.exception))
        s = sample("7/2", "1.5", mutate=lambda s: s["values"][0].pop("abs_err"))
        with self.assertRaises(Refused):
            read_sample(write(self.d, "noerr.json", s))

    def test_masters_from_sample_requires_the_frame(self):
        """review GATE 1 A3: masters_from_sample needs --frame-sha256 for every master; a foreign frame is refused"""
        d0 = "7/2"
        idx = DEFS["m1"]["index"]
        X = write(self.d, "X.json", sample(d0, "1.0", extras=[{"target": idx, "value": "3.25", "abs_err": "1e-60"}]))
        a = write(self.d, "a.json", answer())
        out = Path(self.d) / "m.json"
        rc, _, err = run_tool("masters_from_sample.py", "--answers", a, "--d0", d0, "--out", out, "--master", f"m1={X}")
        self.assertEqual(rc, 2)
        self.assertIn("--frame-sha256", err)
        rc, _, err = run_tool("masters_from_sample.py", "--answers", a, "--d0", d0, "--out", out, "--master", f"m1={X}", "--frame-sha256", "m1=" + "9" * 64)
        self.assertEqual(rc, 2)
        rc, _, err = run_tool("masters_from_sample.py", "--answers", a, "--d0", d0, "--out", out, "--master", f"m1={X}", "--frame-sha256", f"m1={FSHA}")
        self.assertEqual(rc, 0, err)
        m = json.loads(out.read_text())
        from pflowio import content_sha256
        self.assertEqual(m["content_sha256"], content_sha256(m))
        # a sample without a two-precision pass is a production refusal; --measurement records it
        X0 = write(self.d, "X0.json", sample(d0, "1.0", tp=0, extras=[{"target": idx, "value": "3.25", "abs_err": "1e-60"}]))
        rc, _, err = run_tool("masters_from_sample.py", "--answers", a, "--d0", d0, "--out", out, "--master", f"m1={X0}", "--frame-sha256", f"m1={FSHA}")
        self.assertEqual(rc, 2)
        self.assertIn("two-precision", err)
        rc, _, err = run_tool("masters_from_sample.py", "--answers", a, "--d0", d0, "--out", out, "--master", f"m1={X0}", "--frame-sha256", f"m1={FSHA}", "--measurement")
        self.assertEqual(rc, 0, err)
        self.assertTrue(json.loads(out.read_text())["downgrades"])


class Fit(unittest.TestCase):
    """samples of V(eps) = sum_k c_k eps^k (k >= -1) in the E scheme with a known series: the fit recovers it; every contract
    violation refuses the fit"""

    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.d = self.t.name
        mp.dps = 60
        self.c = {-1: mpf(-5), 0: mpf(3) / 7, 1: zeta(3), 2: mpf(2), 3: mpf(-1)}
        self.nodes = ["3999/1000", "4001/1000", "3998/1000", "4002/1000", "3997/1000"]

    def tearDown(self):
        self.t.cleanup()

    def spec(self, files, **kw):
        sp = {"schema": "pflow-fit-spec.v1", "nodes": {n: sha256_file(f) for n, f in zip(self.nodes, files)}, "frame_sha256": FSHA,
              "need_target_digits": 18, "need_abs_err": "1e-30", "target": {"name": "t", "master": [1, 1], "m": 1}, "loops": 2,
              "scheme": "E", "pole_order": 1, "production": True}
        sp.update(kw)
        return sp

    def files(self, bad=None, **kw):
        out = []
        sub = tempfile.mkdtemp(dir=self.d)
        for n in self.nodes:
            e = (4 - rat(n)) / 2
            ef = mpf(e.numerator) / e.denominator
            v = sum(c * ef ** k for k, c in self.c.items())
            s = sample(n, mp.nstr(v, 55), "1e-45", **kw)
            if bad and n == self.nodes[2]:
                s = bad(s)
            out.append(write(sub, n.replace("/", "_") + ".json", s))
        return out

    def test_recovers_and_records(self):
        f = self.files()
        cs, used = fitgrid.fit(self.spec(f), f)
        self.assertLess(abs(cs[1]["c"] - zeta(3)), mpf(10) ** -30)
        self.assertEqual(len(used), 5)
        self.assertTrue(all(len(u["sha256"]) == 64 for u in used))
        self.assertGreater(mpf(cs[1]["abs_err"]), 0)

    def test_contract(self):
        f = self.files()
        sp = self.spec(f)
        for name, files, msg in [("missing", f[:-1], "missing nodes"), ("duplicate", f + [f[0]], "second sample")]:
            with self.assertRaises(Refused) as c:
                fitgrid.fit(sp, files)
            self.assertIn(msg, str(c.exception), name)
        rej = self.files(lambda s: sample(s["d0"], s["values"][0]["value"], "1e-45", accepted=False))
        with self.assertRaises(Refused) as c:
            fitgrid.fit(self.spec(rej), rej)
        self.assertIn("REJECTED", str(c.exception))
        extra = write(self.d, "extra.json", sample("3996/1000", "1.0"))
        with self.assertRaises(Refused):
            fitgrid.fit(sp, f + [extra])
        # pinned samples: another accepted file for the same node (an edit, re-sealed) is not the pinned one
        other = self.files(lambda s: sample(s["d0"], "1.5", "1e-45"))
        with self.assertRaises(Refused) as c:
            fitgrid.fit(sp, other)
        self.assertIn("not the sample pinned", str(c.exception))
        big = self.files(lambda s: sample(s["d0"], s["values"][0]["value"], "1e-10"))
        with self.assertRaises(Refused) as c:
            fitgrid.fit(self.spec(big), big)
        self.assertIn("above the need", str(c.exception))
        for bad, msg in [(dict(sp, frame_sha256=None), "frame_sha256 is required"), (dict(sp, nodes=self.nodes), "pin"),
                         (dict(sp, production=None), "true or false"), (dict(sp, production=0), "true or false"),
                         ({k: v for k, v in sp.items() if k != "production"}, "true or false"),
                         (dict(sp, allow_legacy=True), "replaced by legacy_samples")]:
            with self.assertRaises(Refused) as c:
                fitgrid.fit(bad, f)
            self.assertIn(msg, str(c.exception))

    def test_production_and_legacy(self):
        """review GATE 1 A4 and the two-precision rule: a production contract refuses downgraded samples; legacy_samples only for
        the 12 archived four-loop nodes, never at five loops, never in production"""
        f = self.files(tp=0)
        with self.assertRaises(Refused) as c:
            fitgrid.fit(self.spec(f), f)
        self.assertIn("two-precision", str(c.exception))
        cs, used = fitgrid.fit(self.spec(f, production=False), f)
        self.assertTrue(all(u["downgrades"] for u in used))
        from pflowio import LEGACY_4L_NODES
        arch = next(iter(LEGACY_4L_NODES))
        for kw, msg in [(dict(loops=5, legacy_samples={self.nodes[0]: arch}, production=False), "loops = 5"),
                        (dict(loops=4, legacy_samples={self.nodes[0]: "5" * 64}, production=False), "not one of the 12"),
                        (dict(loops=4, legacy_samples={self.nodes[0]: arch}), "not one of the 12")]:
            with self.assertRaises(Refused) as c:
                fitgrid.fit(self.spec(f, **kw), f)
            self.assertIn(msg, str(c.exception))

    def test_amplification_of_the_production_grid(self):
        """the exact worst-case amplification of c_1 for the 40-node grid grid40-h6 and a fifth-order pole: 10^30.848 (Astra
        review 18 Nov, flow-run/stdout.log)"""
        grid = json.loads((H / "specs" / "grid40-h6.json").read_text())["eps"]
        eps = [rat(x) for x in grid]
        w = fitgrid.cardinal_weights(eps, 1 + 5)
        import math
        amp = sum(abs(x) for x in w)
        self.assertAlmostEqual(math.log10(amp.numerator) - math.log10(amp.denominator), 30.848, places=3)


if __name__ == "__main__":
    unittest.main()
