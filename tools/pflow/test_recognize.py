"""GATE 3 tests: the blind PSLQ protocol on the published controls of the Astra review (18 Nov 2026, recognition.md):
With the CORRECTED search (all inputs rounded to p digits, tolerance n H 10^(1-p), H = 10^12; Astra's fixed 10^(3-p) rejects true
relations of height above ~10^3) the thresholds measured 18 Nov (s1/runs/gate3/scan2-*.json) are:
- GGPPSS int200 finite part after the pole-derived subtraction, nine-element HATTED dictionary: wrong vector at 60 digits, the
  published vector from 70 on (Astra, fixed tolerance: 50 / 60);
- an injected 1e-100 bias is found with the "right" vector at 85 digits but REJECTED by the independent 180-digit confirmation;
- LSS M_{6,2} (rational prefactor removed) eps^3 coefficient (weight <= 10), true 16-element support: wrong at 140, exact from 150;
- its eps^4 coefficient (weight <= 11), 22 elements: wrong at 220 (exact from 230, not tested here: slow).
The TARGETS are built from published formulas with constants that do NOT come from the tool's `mzv`: single zetas from mpmath, the
depth-2 MZVs from the independent polylog quadrature (`mzv_independent`). Only the weight-11 control's zeta(2,1,8) (depth 3) uses
the tool's own evaluation; depth 3 is covered by the stuffle / duality identity tests below.
run: python3 -I -m unittest tools/pflow/test_recognize.py (a few minutes: the quadratures)"""
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mpmath import mp, mpf, zeta, pi  # noqa: E402
import recognize as R  # noqa: E402


def el(name, *terms):
    return {"name": name, "terms": [[str(c), prod] for c, prod in terms]}


def int200():
    """GGPPSS MasterExpansions row 200, finite part, minus the pole-derived -196 zeta(8); hatted zeta(3,5) (GGPPSS order, last index
    outer) = zeta(5,3) outer first - 29/12 zeta(8)"""
    mp.dps = 260
    z = zeta
    z53 = R.mzv_independent((5, 3), 250)
    actual = (mpf(22478) / 23625 * pi ** 8 - mpf(19008) / 5 * z53 - 672 * z(3) ** 2 - 704 * z(3) ** 3 - 7888 * z(3) * z(5)
              + 3108 * z(7) + mpf(1280) / 9 * z(9))
    x = actual + 196 * z(8)
    d = {"schema": R.DICT_SCHEMA, "convention": "MZV outer first", "max_height": 10 ** 12, "maxsteps": 10 ** 6, "elements": [
        el("1", (1, [])), el("z3", (1, [[3]])), el("z5", (1, [[5]])), el("z3^2", (1, [[3], [3]])), el("z7", (1, [[7]])),
        el("z3z5", (1, [[3], [5]])), el("z35hat", (1, [[5, 3]]), (F(-29, 12), [[8]])), el("z9", (1, [[9]])), el("z3^3", (1, [[3], [3], [3]]))]}
    truth = {"z3^2": F(-672), "z7": F(3108), "z3z5": F(-7888), "z35hat": F(-19008, 5), "z9": F(1280, 9), "z3^3": F(-704)}
    return x, d, truth


def lss():
    """LSS M_{6,2} coefficients inside the braces (Astra higher_weight.py; published, arXiv:1108.0732 BC3.tex)"""
    mp.dps = 260
    z = zeta
    n = lambda *ss: R.mzv(tuple(ss), 250)
    names = ["33", "6", "34", "7", "35", "8", "m26", "333", "45", "36", "9", "334", "55", "37", "10", "m37", "335", "56", "47", "38",
             "29", "11", "3m26", "m218"]
    prods = {"33": [[3], [3]], "6": [[6]], "34": [[3], [4]], "7": [[7]], "35": [[3], [5]], "8": [[8]], "m26": [[2, 6]], "333": [[3], [3], [3]],
             "45": [[4], [5]], "36": [[3], [6]], "9": [[9]], "334": [[3], [3], [4]], "55": [[5], [5]], "37": [[3], [7]], "10": [[10]],
             "m37": [[3, 7]], "335": [[3], [3], [5]], "56": [[5], [6]], "47": [[4], [7]], "38": [[3], [8]], "29": [[2], [9]],
             "11": [[11]], "3m26": [[3], [2, 6]], "m218": [[2, 1, 8]]}
    c10 = dict(zip(names[:16], map(F, ["10704", "300", "18138", "156783/2", "214920", "-215670", "59940", "73244/3", "-28884", "13200",
                                       "-1510937/18", "6792", "143960/7", "495583/7", "-1038970/7", "-88260/7"])))
    c11 = dict(zip(names[2:], map(F, ["32112", "146088", "908472", "-395297", "167994", "140784", "-37260", "233700", "166790/3", "109866",
                                     "8380551/28", "26770505/28", "-108421795/56", "-4631493/28", "390194", "-2438995/6", "488786",
                                     "-1567351/2", "3792685", "-100223975/16", "-11366", "205010"])))
    out = {}
    for w, c in ((10, c10), (11, c11)):
        d = {"schema": R.DICT_SCHEMA, "convention": "MZV outer first", "max_height": 10 ** 12, "maxsteps": 10 ** 6,
             "elements": [el(k, (1, prods[k])) for k in c]}
        # target from INDEPENDENT constants (depth 3 only: the tool's own evaluation, see the module docstring)
        def ind(ss):
            v = R.mzv_independent(tuple(ss), 250)
            return v if v is not None else R.mzv(tuple(ss), 250)
        mp.dps = 260
        vals = []
        for k in c:
            t = mpf(1)
            for ss in prods[k]:
                t *= ind(ss)
            vals.append(t)
        x = sum(mpf(v.numerator) / v.denominator * b for v, b in zip(c.values(), vals))
        out[w] = (x, d, c)
    return out


class Constants(unittest.TestCase):
    def test_identities_and_the_independent_evaluation(self):
        """stuffle, duality and pi^6/7! identities at 200 digits (they check the word convention and the algebra, which the two
        split points of one algorithm cannot), and the independent polylog quadrature for depth-2 values"""
        mp.dps = 220
        z = lambda *ss: R.mzv(tuple(ss), 210)
        eps = mpf(10) ** -200
        self.assertLess(abs(z(3, 3) - (zeta(3) ** 2 - zeta(6)) / 2), eps)
        self.assertLess(abs(z(4, 4) - (zeta(4) ** 2 - zeta(8)) / 2), eps)
        self.assertLess(abs(z(5, 3) + z(3, 5) + zeta(8) - zeta(3) * zeta(5)), eps)
        self.assertLess(abs(z(7, 3) + z(3, 7) + zeta(10) - zeta(3) * zeta(7)), eps)
        self.assertLess(abs(z(2, 1, 1) - zeta(4)), eps)
        self.assertLess(abs(z(2, 2, 2) - pi ** 6 / 5040), eps)
        self.assertLess(abs(z(3, 1, 2) - z(2, 3, 1)), eps)  # duality
        for ss in [(5, 3), (7, 3), (2, 6)]:
            self.assertLess(abs(R.mzv(ss, 150) - R.mzv_independent(ss, 150)), mpf(10) ** -145, ss)


    def test_two_methods_and_known_values(self):
        mp.dps = 160
        self.assertLess(abs(R.mzv((2, 1), 150) - zeta(3)), mpf(10) ** -148)
        self.assertLess(abs(R.mzv((5, 3), 150) - mpf("0.037707672984847544011304782")), mpf(10) ** -26)
        for ss in [(3, 5), (2, 1, 8), (3, 7)]:
            self.assertLess(abs(R.mzv(ss, 150, (1, 2)) - R.mzv(ss, 150, (1, 3))), mpf(10) ** -148)
        self.assertEqual([len(R.hoffman_dictionary(W)["elements"]) for W in (9, 10, 11)], [20, 27, 36])


class Int200(unittest.TestCase):
    def test_scan_and_protocol(self):
        x, d, truth = int200()
        rows = R.scan(x, d, [60, 70], truth=truth, const_dps=200)
        self.assertEqual([r["status"] for r in rows], ["other", "exact"])
        # protocol: search at 70 on the clean value, confirmation on an independent 180-digit value
        ok = R.recognize(x, x, mpf(10) ** -178 * abs(x), d, 70, const_dps=200)
        self.assertEqual(ok["verdict"], "ACCEPTED", ok)
        # a 1e-100 bias: the search at 85 still returns the published vector; the independent 180-digit value of the SAME biased
        # quantity rejects it (residual 1e-100 >> its uncertainty)
        xb = x + mpf(10) ** -100
        rel, _ = R.search(xb, R.constants(d, 200), 85, d["max_height"], d["maxsteps"])
        self.assertEqual({d["elements"][i]["name"]: c for i, c in rel.items()}, truth)
        bad = R.recognize(xb, xb, mpf(10) ** -178 * abs(x), d, 85, const_dps=200)
        self.assertEqual(bad["verdict"], "REJECTED", bad)
        self.assertIn("residual", bad["confirmation"]["reason"])
        # too few independent digits: 120-digit confirmation of a 85-digit search is refused (excess 35 < 40)
        few = R.recognize(x, x, mpf(10) ** -120 * abs(x), d, 85, const_dps=200)
        self.assertEqual(few["verdict"], "REJECTED")
        self.assertIn("excess", few["confirmation"]["reason"])


class HigherWeight(unittest.TestCase):
    def test_weight_10_and_11(self):
        L = lss()
        x, d, c = L[10]
        rows = R.scan(x, d, [140, 150], truth=c, const_dps=220)
        self.assertEqual([r["status"] for r in rows], ["other", "exact"])
        x, d, c = L[11]
        rows = R.scan(x, d, [220], truth=c, const_dps=260)
        self.assertNotEqual(rows[0]["status"], "exact")


class Cli(unittest.TestCase):
    def test_command_line_on_the_int200_control(self):
        import json, subprocess, tempfile
        x, d, truth = int200()
        with tempfile.TemporaryDirectory() as t:
            dp = Path(t) / "d.json"
            dp.write_text(json.dumps(d))
            mp.dps = 260
            run = lambda *a: subprocess.run([sys.executable, "-I", "-B", str(Path(__file__).with_name("recognize.py")), *a], capture_output=True, text=True)
            r = run("--dict", str(dp), "--dict-sha256", R.dict_sha(d), "--p1", "70", "--x-a", mp.nstr(x, 120), "--x-b", mp.nstr(x, 200), "--err-b", mp.nstr(abs(x) * mpf(10) ** -178, 3))
            self.assertTrue(r.stdout.startswith("ACCEPTED"), r.stdout + r.stderr)
            r = run("--dict", str(dp), "--dict-sha256", R.dict_sha(d), "--p1", "70", "--x-a", mp.nstr(x, 120), "--x-b", mp.nstr(x + mpf(10) ** -100, 200), "--err-b", mp.nstr(abs(x) * mpf(10) ** -178, 3))
            self.assertTrue(r.stdout.startswith("REJECTED"), r.stdout + r.stderr)


class CliFits(unittest.TestCase):
    """the command line with fit files (review G1, G2): sealed fits of the same coefficient on disjoint nodes, fit A with enough
    digits, a pinned dictionary"""

    def test_fit_files_and_dictionary_pin(self):
        import json, subprocess, tempfile, copy
        from pflowio import seal
        x, d, truth = int200()
        mp.dps = 260
        spec = {"schema": "pflow-fit-spec.v1", "target": {"name": "t", "master": [1]}, "frame_sha256": "f" * 64, "scheme": "G", "loops": 5,
                "pole_order": 5, "production": True}
        def fit(nodes, val, err, samples):
            return {"schema": "pflow-fit.v1", "spec": dict(spec, nodes={n: "0" * 64 for n in nodes}),
                    "used_samples": [{"sha256": h} for h in samples], "coefficients": {"1": {"value": mp.nstr(val, 230), "abs_err": err}}}
        A = fit(["1/2", "1/3"], x, mp.nstr(abs(x) * mpf(10) ** -110, 3), ["a1", "a2"])
        B = fit(["1/4", "1/5"], x, mp.nstr(abs(x) * mpf(10) ** -178, 3), ["b1", "b2"])
        with tempfile.TemporaryDirectory() as t:
            T = Path(t)
            dp = T / "d.json"
            dp.write_text(json.dumps(d))
            sha = R.dict_sha(d)
            def run(a, b, *extra, pin=sha):
                (T / "a.json").write_text(json.dumps(a))
                (T / "b.json").write_text(json.dumps(b))
                args = ["--dict", str(dp), "--p1", "70", "--fit-a", str(T / "a.json"), "--fit-b", str(T / "b.json"), "--coef", "1", *extra]
                if pin:
                    args += ["--dict-sha256", pin]
                r = subprocess.run([sys.executable, "-I", "-B", str(Path(__file__).with_name("recognize.py")), *args], capture_output=True, text=True)
                return r.stdout + r.stderr
            self.assertTrue(run(seal(copy.deepcopy(A)), seal(copy.deepcopy(B))).startswith("ACCEPTED"))
            ed = seal(copy.deepcopy(B))
            ed["coefficients"]["1"]["abs_err"] = "1e-300"
            for a, b, msg, kw in [
                (A, B, "unsealed", {}),
                (seal(copy.deepcopy(A)), ed, "edited", {}),
                (seal(copy.deepcopy(A)), seal(dict(copy.deepcopy(B), spec=dict(B["spec"], target={"name": "u", "master": [2]}))), "differ in target", {}),
                (seal(copy.deepcopy(A)), seal(dict(copy.deepcopy(B), spec=dict(B["spec"], nodes={"1/2": "0" * 64}))), "share nodes", {}),
                (seal(fit(["1/2"], x, mp.nstr(abs(x) * mpf(10) ** -50, 3), ["a1"])), seal(copy.deepcopy(B)), "informative digits", {}),
                (seal(copy.deepcopy(A)), seal(copy.deepcopy(B)), "not the pinned", {"pin": "9" * 64}),
                (seal(copy.deepcopy(A)), seal(copy.deepcopy(B)), "needs --dict-sha256", {"pin": None}),
            ]:
                out = run(a, b, **kw)
                self.assertIn("REFUSED", out, msg)
                self.assertIn(msg, out)


if __name__ == "__main__":
    unittest.main()
