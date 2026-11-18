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
from pflowio import Refused, eval_exact, read_masters, read_sample, rat  # noqa: E402
import fitgrid  # noqa: E402


def g11(eps):
    return gamma(1 + eps) * gamma(1 - eps) ** 2 / gamma(2 - 2 * eps)


def run_tool(name, *args):
    p = subprocess.run([sys.executable, "-I", "-B", str(H / name), *map(str, args)], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


class Exact(unittest.TestCase):
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

    def masters(self, d0, err="1e-50"):
        from pflowio import digest
        return write(self.d, f"m-{d0.replace('/', '_')}.json", {"schema": "pflow-masters.v1", "d0": d0, "master_definitions_sha256": digest(DEFS),
                                                              "values": {"m1": {"value": "2.5", "abs_err": err}}, "sources": []})

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
        mf = write(self.d, "mbad.json", m)
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


def sample(d0, value, err="1e-40", accepted=True, frame="f" * 64, target=30.0, **kw):
    s = {"schema": "pflow-route-b.v2", "accepted": accepted, "acceptance": {"accepted": accepted, "failures": [] if accepted else ["x"]},
         "accept_diagnostics": False, "d0": d0, "frame_sha256": frame, "target_digits": target,
         "values": [{"m": 1, "master": [1, 1], "value": value, "abs_err": err}], "extras": []}
    s.update(kw)
    return s


class Samples(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.d = self.t.name

    def tearDown(self):
        self.t.cleanup()

    def test_refusals(self):
        ok = write(self.d, "ok.json", sample("7/2", "1.5"))
        self.assertEqual(read_sample(ok, d0="7/2", frame_sha="f" * 64, need_target=20)["values"][(1, (1, 1))][0], "1.5")
        for name, s, kw, msg in [
            ("rej", sample("7/2", "1.5", accepted=False), {}, "REJECTED"),
            ("nov", {k: v for k, v in sample("7/2", "1.5").items() if k != "accepted"}, {}, "no acceptance verdict"),
            ("diag", sample("7/2", "1.5", accept_diagnostics=True), {}, "accept-diagnostics"),
            ("d0", sample("15/4", "1.5"), {"d0": "7/2"}, "required d0"),
            ("frame", sample("7/2", "1.5", frame="e" * 64), {"frame_sha": "f" * 64}, "frame_sha256"),
            ("tgt", sample("7/2", "1.5", target=12.0), {"need_target": 18}, "needed"),
            ("v1", {"schema": "pflow-route-b.v1", "accepted": True, "d0": "7/2", "values": [], "config": {"w": 40}}, {}, "schema"),
        ]:
            with self.assertRaises(Refused) as c:
                read_sample(write(self.d, name + ".json", s), **kw)
            self.assertIn(msg, str(c.exception), name)
        s = sample("7/2", "1.5")
        del s["values"][0]["abs_err"]
        with self.assertRaises(Refused):
            read_sample(write(self.d, "noerr.json", s))


class Fit(unittest.TestCase):
    """samples of V(eps) = sum_k c_k eps^k (k >= -1) in the E scheme with a known series: the fit recovers it; every contract
    violation refuses the fit"""

    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.d = self.t.name
        mp.dps = 60
        self.c = {-1: mpf(-5), 0: mpf(3) / 7, 1: zeta(3), 2: mpf(2), 3: mpf(-1)}
        self.nodes = ["3999/1000", "4001/1000", "3998/1000", "4002/1000", "3997/1000"]
        self.spec = {"schema": "pflow-fit-spec.v1", "nodes": self.nodes, "frame_sha256": "f" * 64, "need_target_digits": 18,
                     "need_abs_err": "1e-30", "target": {"name": "t", "master": [1, 1], "m": 1}, "loops": 2, "scheme": "E", "pole_order": 1}

    def tearDown(self):
        self.t.cleanup()

    def files(self, bad=None):
        out = []
        sub = tempfile.mkdtemp(dir=self.d)
        for n in self.nodes:
            e = (4 - rat(n)) / 2
            ef = mpf(e.numerator) / e.denominator
            v = sum(c * ef ** k for k, c in self.c.items())
            s = sample(n, mp.nstr(v, 55), "1e-45")
            if bad and n == self.nodes[2]:
                s = bad(s)
            out.append(write(sub, n.replace("/", "_") + ".json", s))
        return out

    def test_recovers_and_records(self):
        cs, used = fitgrid.fit(self.spec, self.files())
        self.assertLess(abs(cs[1]["c"] - zeta(3)), mpf(10) ** -30)
        self.assertEqual(len(used), 5)
        self.assertTrue(all(len(u["sha256"]) == 64 for u in used))
        self.assertGreater(mpf(cs[1]["abs_err"]), 0)

    def test_contract(self):
        f = self.files()
        for name, files, msg in [("missing", f[:-1], "missing nodes"), ("duplicate", f + [f[0]], "second sample"),
                                 ("rejected", self.files(lambda s: dict(s, accepted=False)), "REJECTED")]:
            with self.assertRaises(Refused) as c:
                fitgrid.fit(self.spec, files)
            self.assertIn(msg, str(c.exception), name)
        extra = write(self.d, "extra.json", sample("3996/1000", "1.0"))
        with self.assertRaises(Refused):
            fitgrid.fit(self.spec, f + [extra])

        def big_err(s):
            s["values"][0]["abs_err"] = "1e-10"
            return s
        with self.assertRaises(Refused) as c:
            fitgrid.fit(self.spec, self.files(big_err))
        self.assertIn("above the need", str(c.exception))
        sp = dict(self.spec, frame_sha256=None)
        with self.assertRaises(Refused):
            fitgrid.fit(sp, f)

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
