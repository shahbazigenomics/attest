"""Confidence intervals and tests for the validation figures reported in the paper.

Run from the attest project root:  python validation/statistics/validation_statistics.py
Needs scipy. Reads only CSVs already produced by the validation scripts; writes
validation/statistics/validation_statistics.csv and prints a markdown table.
Proportions: Wilson 95% CI, plus the exact (Clopper-Pearson) bound where a count
is 0 or n, since Wilson is not the conventional choice at the boundary.
"""
import csv, collections, math
from scipy.stats import beta, chi2_contingency, norm

Z = norm.ppf(0.975)

def wilson(k, n):
    p = k / n; d = 1 + Z**2 / n
    c = (p + Z**2 / (2 * n)) / d; h = Z * math.sqrt(p * (1 - p) / n + Z**2 / (4 * n * n)) / d
    return max(0, c - h), min(1, c + h)

def exact(k, n):
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return lo, hi

rows = []
def prop(group, what, k, n, source):
    lo, hi = wilson(k, n); elo, ehi = exact(k, n)
    rows.append(dict(group=group, measure=what, k=k, n=n, estimate=round(100 * k / n, 1),
                     wilson_lo=round(100 * lo, 1), wilson_hi=round(100 * hi, 1),
                     exact_lo=round(100 * elo, 1), exact_hi=round(100 * ehi, 1), source=source))

V = "validation/"
# attest: value scale on full real matrices (13 transforms x 3 species)
r = list(csv.DictReader(open(V + "real_data_results.csv")))
prop("attest", "value-scale verdicts correct, 3 species x 13 transforms", sum(x["agree"] == "TRUE" for x in r), len(r), "real_data_results.csv")

# attest: depth sweep
r = list(csv.DictReader(open(V + "depth_sweep.csv")))
for d in sorted({float(x["depth"]) for x in r}, reverse=True):
    s = [x for x in r if float(x["depth"]) == d]
    prop("attest", f"depth sweep: verdicts correct at {d:,.0f} reads/sample", sum(x["correct"] == "TRUE" for x in s), len(s), "depth_sweep.csv")
s = [x for x in r if float(x["depth"]) >= 1e6]
prop("attest", "depth sweep: verdicts correct at >= 1M reads/sample", sum(x["correct"] == "TRUE" for x in s), len(s), "depth_sweep.csv")
s = [x for x in r if x["expected"] == "PERMITTED"]
prop("attest", "depth sweep: false alarms on raw counts, all depths", sum(x["got"] != "PERMITTED" for x in s), len(s), "depth_sweep.csv")

# attest: simulated raw matrices
r = list(csv.DictReader(open(V + "simulation_specificity.csv")))
n = sum(int(x["PERMITTED"]) + int(x["CAUTION"]) + int(x["NOT"]) + int(x["UNKNOWN"]) for x in r)
prop("attest", "false alarms on simulated raw matrices (9 regimes)", n - sum(int(x["PERMITTED"]) for x in r), n, "simulation_specificity.csv")

# attest: GEO study (original run, 341 audited files)
g = [x for x in csv.DictReader(open(V + "geo_study/results/files.csv")) if x["status"] == "audited"]
nonraw = [x for x in g if x["truth_value"] in ("depth removed", "log-transformed")]
raw = [x for x in g if x["truth_value"] == "raw counts"]
prop("attest", "GEO: not-raw files not called PERMITTED (sensitivity, lenient)", sum(x["v_value"] != "PERMITTED" for x in nonraw), len(nonraw), "geo_study/results/files.csv")
prop("attest", "GEO: not-raw files called NOT PERMITTED (sensitivity, strict)", sum(x["v_value"] == "NOT PERMITTED" for x in nonraw), len(nonraw), "geo_study/results/files.csv")
prop("attest", "GEO: raw-count files called PERMITTED (specificity)", sum(x["v_value"] == "PERMITTED" for x in raw), len(raw), "geo_study/results/files.csv")
prop("attest", "GEO: files called NOT PERMITTED (headline)", sum(x["v_value"] == "NOT PERMITTED" for x in g), len(g), "geo_study/results/files.csv")
DET = {"raw counts", "depth removed", "log-transformed", "count scale, not whole numbers", "estimated counts"}
det = [x for x in g if x["truth_value"] in DET]; und = [x for x in g if x["truth_value"] not in DET]
prop("attest", "GEO: files with an NCBI-based truth (truth coverage)", len(det), len(g), "geo_study/results/files.csv")

# selection bias: verdict distribution, truth determined vs not
levels = ["PERMITTED", "CAUTION", "NOT PERMITTED", "UNKNOWN"]
tab = [[sum(x["v_value"] == v for x in det) for v in levels], [sum(x["v_value"] == v for x in und) for v in levels]]
chi = chi2_contingency(tab)
tab3 = [t[:3] for t in tab]          # without UNKNOWN, which cannot have a truth by construction
chi3 = chi2_contingency(tab3)

# run-to-run agreement (frozen-code re-run)
c = [x for x in csv.DictReader(open(V + "geo_study/results_rerun_2026-11/concordance_vs_2026-09-22.csv"))
     if x["status_2026-09-22"] == "audited" and x["status_rerun"] == "audited"]
a = [x["v_value_2026-09-22"] for x in c]; b = [x["v_value_rerun"] for x in c]
po = sum(i == j for i, j in zip(a, b)) / len(a)
cats = sorted(set(a) | set(b))
pe = sum((a.count(k) / len(a)) * (b.count(k) / len(b)) for k in cats)
kappa = (po - pe) / (1 - pe)
prop("attest", "GEO re-run: value-scale verdict identical", sum(i == j for i, j in zip(a, b)), len(a), "concordance_vs_2026-09-22.csv")

# admissible, CEPH 1463 (v0.1.2; tests/test_public_ceph.py and `admissible identity --json`)
prop("admissible", "CEPH: sex calls wrong among samples called", 0, 16, "admissible identity, CEPH 1463")
prop("admissible", "CEPH: unrelated pairs called related", 0, 11, "admissible identity, CEPH 1463")
prop("admissible", "CEPH: first-degree pairs recognised as related", 81, 81, "admissible identity, CEPH 1463")
prop("admissible", "CEPH: false duplicate pairs", 0, 136, "admissible identity, CEPH 1463")
prop("admissible", "ANNOVAR reader: genotypes reconstructed correctly", 120, 120, "tests/test_annovar_real.py")

out = "validation/statistics/validation_statistics.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print("| tool | measure | k/n | % | Wilson 95% CI | exact 95% CI |\n|---|---|---|---|---|---|")
for x in rows:
    print(f"| {x['group']} | {x['measure']} | {x['k']}/{x['n']} | {x['estimate']} | {x['wilson_lo']}-{x['wilson_hi']} | {x['exact_lo']}-{x['exact_hi']} |")
print("\nSelection bias, attest verdict by truth determined vs not (rows: determined, undetermined; cols:", levels, ")")
print(tab, f"chi2={chi[0]:.2f}, df={chi[2]}, p={chi[1]:.2g}")
print("excluding UNKNOWN:", tab3, f"chi2={chi3[0]:.2f}, df={chi3[2]}, p={chi3[1]:.2g}")
print(f"\nRun-to-run value-scale agreement: {po*100:.1f}% (n={len(a)}), Cohen's kappa = {kappa:.3f}")
