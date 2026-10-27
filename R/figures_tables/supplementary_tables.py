# =============================================================================
# supplementary_tables.py
# Builds Tables S1 to S10 and S12 of the paper from the result files of the
# scripts in R/validation/, in the order of first citation. Table S11 (mapped
# area by distance class) comes from the two CSV exports of
# gee/11_distance_to_field_samples.js.
#
# Usage (from the repository root):
#   python3 R/figures_tables/supplementary_tables.py result result/supplement
# Output: one CSV per table in result/supplement/tables_csv/
# =============================================================================
import glob, json, os, sys
import pandas as pd

res_dir = sys.argv[1] if len(sys.argv) > 1 else "result"
out_dir = sys.argv[2] if len(sys.argv) > 2 else "result/supplement"
os.makedirs(os.path.join(out_dir, "tables_csv"), exist_ok=True)

def f(name):
    p = sorted(glob.glob(os.path.join(res_dir, "*" + name)))
    if not p:
        raise FileNotFoundError(name)
    return pd.read_csv(p[0])

M = "−"            # minus sign, as in the main text
def num(x, d=2):
    if pd.isna(x):
        return "–"
    s = f"{x:.{d}f}"
    if s.startswith("-"):
        s = M + s[1:]
    if s in (M + "0." + "0" * d, M + "0"):
        s = s[1:]
    return s
def rng(lo, hi, d=2):
    """Interval text: (a–b) when both are positive, (a to b) when a sign is negative."""
    a, b = num(lo, d), num(hi, d)
    return f"({a} to {b})" if (M in a or M in b) else f"({a}–{b})"
def mr(m, lo, hi, d=2):
    return f"{num(m, d)} {rng(lo, hi, d)}"
def pct(x, d=1):
    s = num(x, d)
    return ("+" + s) if (not s.startswith(M) and float(x) > 0.05) else s
def intc(x):
    return f"{int(round(x)):,}"

REGIONS = ["Ayau", "Bintan", "Karimunjawa", "Komodo", "Manado", "Rote"]
GSE, PX, BL = "GSE bands only", "Proxy workflow", "Mean-value baseline"
tables = []

def add(tid, caption, header, rows, widths, note=None, landscape=False,
        bold=None, highlight=None, csv_rows=None):
    tables.append(dict(id=tid, caption=caption, header=header, rows=rows, widths=widths,
                       note=note, landscape=landscape, bold=bold or [], highlight=highlight or []))
    path = os.path.join(out_dir, "tables_csv", tid.replace(" ", "") + ".csv")
    if csv_rows:                      # first row of csv_rows is the header
        pd.DataFrame(csv_rows).to_csv(path, index=False, header=False, encoding="utf-8-sig")
    else:
        pd.DataFrame(rows, columns=header[-1]).to_csv(path, index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------------------
# Table S1  Samples and accuracy of the proxy models (random vs pixel-grouped)
# ---------------------------------------------------------------------------
fp = f("proxy_stage_MC_summary_fullPipeline.csv")
pt = f("proxy_stage_MC_summary_pointTarget.csv")
def st(df, stage, metric, mode, d=2):
    r = df[(df.Stage == stage) & (df.Metric == metric) & (df.Mode == mode)].iloc[0]
    return mr(r["mean"], r["lower"], r["upper"], d)
s1 = [
    ("Seagrass presence", "5,408 (3,876 present, 1,532 absent)", "Overall accuracy", fp, "PA", "Accuracy", 2),
    ("", "", "F1-score (presence)", fp, "PA", "F1", 2),
    ("Leaf morphology", "3,599 (3,119 short leaf, 263 mixed leaf, 217 long leaf)", "Overall accuracy", fp, "MORPH3", "Accuracy", 2),
    ("", "", "Kappa", fp, "MORPH3", "Kappa", 2),
    ("Percent cover", "5,395 (3,831 with seagrass)", "R²", pt, "SPC_presence", "R2", 2),
    ("", "", "RMSE (% cover)", pt, "SPC_presence", "RMSE", 1),
    ("Above-ground biomass", "3,676", "R²", pt, "AGB", "R2", 2),
    ("", "", "RMSE (gDW m⁻²)", pt, "AGB", "RMSE", 1),
    ("Carbon index", "3,599", "R²", fp, "CINDEX", "R2", 2),
    ("", "", "RMSE", fp, "CINDEX", "RMSE", 2),
]
rows = [[a, b, c, st(df, s, m, "Random rows", d), st(df, s, m, "Pixel-grouped", d)]
        for a, b, c, df, s, m, d in s1]
add("Table S1",
    "Samples and accuracy of the five proxy models. Accuracy is the mean of ten Monte Carlo repeats "
    "with 30% of the samples held out, with the 2.5th and 97.5th percentiles in parentheses. The "
    "random-sample split held out single samples. The pixel-grouped split held out all samples of a "
    "pixel together. Percent cover was scored only on samples with seagrass present.",
    [["Proxy model", "Samples (n)", "Metric", "Random-sample split", "Pixel-grouped split"]],
    rows, [16, 30, 18, 18, 18],
    note="Long leaf, Enhalus acoroides only; short leaf, short-leaf species only; mixed leaf, "
         "short-leaf species with E. acoroides. RMSE, root mean square error.")

# ---------------------------------------------------------------------------
# Table S2  Species weights of the carbon index
# ---------------------------------------------------------------------------
# Values from the shoot carbon spreadsheet ("last update" block): AGBshoot_Mean,
# %OCC from Rahmawati et al., and AGCshoot_Mean = the weight in R/workflow/06_carbon_index_field.R
w = [("Enhalus acoroides", "Ea", 2.369, 35, 0.829), ("Thalassia hemprichii", "Th", 0.973, 33, 0.321),
     ("Cymodocea rotundata", "Cr", 0.400, 38, 0.152), ("Cymodocea serrulata", "Cs", 0.592, 32, 0.190),
     ("Syringodium isoetifolium", "Si", 0.171, 31, 0.053), ("Halodule uninervis", "Hu", 0.044, 32, 0.014),
     ("Halophila ovalis", "Ho", 0.064, 27, 0.017), ("Halodule pinifolia", "Hp", 0.076, 40, 0.031),
     ("Thalassodendron ciliatum", "Tc", 0.277, 39, 0.108), ("Halophila minor", "Hm", None, 24, 0.0),
     ("Halophila sulawesii", "Hs", None, None, 0.0), ("Halophila decipiens", "Hd", None, None, 0.0)]
rows = [[sp, c, "–" if dw is None else f"{dw:.3f}", "–" if oc is None else str(oc), f"{v:.3f}"]
        for sp, c, dw, oc, v in w]
hl = []
add("Table S2",
    "Species weights used in the carbon index. The weight is the mean above-ground carbon per shoot, "
    "calculated as the mean dry weight of one shoot multiplied by the organic carbon content (OCC) of "
    "the species. Shoot dry weights are means from five Indonesian field surveys (Rote 2021, "
    "Karimunjawa 2011 and 2018, Bintan 2019 and Selayar 2015). OCC values are from Rahmawati, Hernawan "
    "and Rustam (2019).",
    [["Species", "Code", "Mean shoot dry weight (gDW shoot⁻¹)", "OCC (%)",
      "Weight (gC shoot⁻¹)"]],
    rows, [32, 10, 22, 14, 22], highlight=hl,
    note="Each shoot dry weight is the dry weight of a single shoot. Weights were calculated before "
         "rounding. Hm, Hs and Hd had no shoot data and were given a weight of 0.")

# ---------------------------------------------------------------------------
# Table S3  Random Forest settings
# ---------------------------------------------------------------------------
rows = [
    ["Seagrass presence", "Presence (0 or 1)", "GSE + Env", "700", "9", "3", "0.8"],
    ["Leaf morphology", "Three leaf classes", "GSE + Env", "300", "6", "7", "0.7"],
    ["Percent cover", "Total cover (%)", "GSE + Env + 3 morphology probabilities", "300", "10", "9", "0.7"],
    ["Above-ground biomass", "AGB point estimate (gDW m⁻²)", "GSE + Env", "300", "8", "9", "0.6"],
    ["Carbon index", "Carbon index (0–1)", "GSE + Env", "500", "6", "9", "0.6"],
    ["AGC, proxy workflow", "AGC point estimate (gC m⁻²)", "GSE + Env + 7 proxy layers", "900", "5", "9", "0.6"],
    ["AGC, GSE-only model", "AGC point estimate (gC m⁻²)", "GSE + Env", "900", "5", "9", "0.6"],
]
add("Table S3",
    "Random Forest settings of the five proxy models and the AGC models. Settings were chosen by a grid "
    "search over five random 70/30 splits and were kept fixed in all validation designs.",
    [["Model", "Response", "Predictors", "Trees", "mtry", "Min. node size", "Sample fraction"]],
    rows, [20, 22, 24, 8, 7, 9, 10],
    note="GSE, the 64 bands of the Google Satellite Embeddings; Env, water depth and distance to land. "
         "The 7 proxy layers are the presence probability, the three morphology class probabilities, "
         "percent cover, AGB and the carbon index.")

# ---------------------------------------------------------------------------
# Table S4  Proxy models refitted vs trained on all samples
# ---------------------------------------------------------------------------
it = f("agc_proxy_sensitivity_iterations_ext.csv")
rows = []
for split, lab in [("Pixel-grouped", "Pixel hold-out"), ("Random rows", "Random samples")]:
    for prox, plab in [("Refitted", "Refitted with training samples only"),
                       ("Trained on all samples", "Trained on all samples")]:
        x = it[(it.Split == split) & (it.Proxies == prox)]
        q = lambda c, d: mr(x[c].mean(), x[c].quantile(.025), x[c].quantile(.975), d)
        rows.append([lab if prox == "Refitted" else "", plab, q("RMSE", 1), q("R2", 2), q("MEC", 2)])
add("Table S4",
    "AGC accuracy of the proxy workflow when the proxy models were refitted with the training samples "
    "only, and when they were trained on all samples, including the test samples. Values are means of "
    "ten repeats, with the 2.5th and 97.5th percentiles in parentheses.",
    [["Split", "Proxy models", "RMSE (gC m⁻²)", "R²", "MEC"]],
    rows, [18, 31, 17, 17, 17],
    note="Only the refitted runs are leak-free. The runs with proxy models trained on all samples show "
         "how much the accuracy is overstated when this step is skipped.")

# ---------------------------------------------------------------------------
# Table S5  Pixel sharing
# ---------------------------------------------------------------------------
px = f("agc_pixel_sharing_fullPipeline.csv").set_index("Region")
rows = [[r, intc(px.loc[r, "n_quadrats"]), intc(px.loc[r, "n_pixels"]),
         intc(px.loc[r, "max_quadrats_per_pixel"]), num(px.loc[r, "pct_quadrats_shared"], 1)]
        for r in REGIONS + ["All regions"]]
add("Table S5",
    "AGC quadrats that share a 10-m pixel with another quadrat. A pixel is defined by identical GSE "
    "values in the same year.",
    [["Region", "Quadrats", "Pixels", "Max. quadrats per pixel", "Quadrats sharing a pixel (%)"]],
    rows, [24, 16, 16, 20, 24], bold=[(6, c) for c in range(5)])

# ---------------------------------------------------------------------------
# Table S6  Accuracy per region, pixel hold-out
# ---------------------------------------------------------------------------
pi = f("agc_nestedMC_iterations_pointTarget.csv")
rows = []
for reg in REGIONS + ["All regions"]:
    for k, (mod, mlab) in enumerate([(PX, "Proxy workflow"), (GSE, "GSE-only"), (BL, "Baseline")]):
        x = pi[(pi.Scope == reg) & (pi.Model == mod)]
        r2 = "–" if mod == BL else f"{num(x.R2.mean())} {rng(x.R2.min(), x.R2.max())}"
        rows.append([reg if k == 0 else "", intc(x.n_test.mean()) if k == 0 else "", mlab,
                     num(x.RMSE.mean(), 1), r2,
                     mr(x.MEC.mean(), x.MEC.quantile(.025), x.MEC.quantile(.975)),
                     pct(x.Bias_pct.mean())])
add("Table S6",
    "AGC accuracy in each region in the pixel hold-out design. Values are means of ten repeats. "
    "Parentheses give the lowest and highest R² and the 2.5th and 97.5th percentiles of MEC across "
    "repeats. n is the mean number of test samples per repeat. Bias is the bias of the mean AGC of "
    "the test samples.",
    [["Region", "n", "Model", "RMSE (gC m⁻²)", "R² (min to max)", "MEC", "Bias (%)"]],
    rows, [16, 8, 17, 13, 18, 18, 10])

# ---------------------------------------------------------------------------
# Table S7  Proxy values in training and test samples, and proxy importance
# ---------------------------------------------------------------------------
cr = f("agc_crossfit_proxy_r2_ext.csv")
cs = f("agc_crossfit_summary_ext.csv")
names = [("PA_prob", "Presence probability"), ("tSPC_pred", "Percent cover"),
         ("tAGB_pred", "AGB"), ("carbon_index", "Carbon index")]
sets = ["Training, in-sample proxies", "Training, cross-fitted proxies", "Test"]
rows = []
for code, lab in names:
    row = [lab]
    for des in ["Pixel-grouped", "Block 250 m"]:
        for s in sets:
            v = cr[(cr.Design == des) & (cr.Proxy == code) & (cr.Set == s)].R2_with_AGC.iloc[0]
            row.append(num(v))
    rows.append(row)
imp = ["Share of proxy layers in importance (%)"]
for des in ["Pixel-grouped", "Block 250 m"]:
    a = cs[(cs.Design == des) & (cs.Scope == "All regions") & (cs.Model == PX)].share_proxy_imp.iloc[0]
    b = cs[(cs.Design == des) & (cs.Scope == "All regions") &
           (cs.Model == "Proxy workflow, cross-fitted")].share_proxy_imp.iloc[0]
    imp += [num(100 * a, 0), num(100 * b, 0), "–"]
rows.append(imp)
hdr = [[{"text": "Proxy layer", "rowSpan": 2}, {"text": "Pixel hold-out", "colSpan": 3},
        {"text": "250 m blocks", "colSpan": 3}],
       ["Training, in-sample", "Training, cross-fitted", "Test",
        "Training, in-sample", "Training, cross-fitted", "Test"]]
add("Table S7",
    "Share of the AGC variation explained by each proxy layer (R²) in the training and test samples, "
    "and share of the proxy layers in the predictor importance of the AGC model. In-sample proxy values "
    "come from proxy models fitted to the same training samples. Cross-fitted values come from proxy "
    "models fitted without the group of the sample. Values are means of ten repeats.",
    hdr, rows, [28, 12, 12, 12, 12, 12, 12],
    note="Importance is the increase in node purity. The morphology class probabilities are part of the "
         "proxy share but are not listed as rows.",
    csv_rows=[["Proxy layer", "Pixel in-sample", "Pixel cross-fitted", "Pixel test",
               "Block in-sample", "Block cross-fitted", "Block test"]] + rows)

# ---------------------------------------------------------------------------
# Table S8  Block sizes
# ---------------------------------------------------------------------------
bk = f("agc_block_summary_ext.csv")
designs = [("Pixel-grouped (PART D)", "Pixel hold-out"), ("Block 100 m", "100 m blocks"),
           ("Block 250 m", "250 m blocks"), ("Block 500 m", "500 m blocks"),
           ("Block 1000 m", "1000 m blocks")]
cs_px = cs[(cs.Design == "Pixel-grouped") & (cs.Scope == "All regions")].set_index("Model")
rows, bold = [], []
for i, (d, lab) in enumerate(designs):
    a = bk[(bk.Design == d) & (bk.Scope == "All regions")].set_index("Model")
    if d.startswith("Pixel"):
        mg = mr(a.loc[GSE, "MEC_mean"], cs_px.loc[GSE, "MEC_lower"], cs_px.loc[GSE, "MEC_upper"])
        mb = mr(a.loc[BL, "MEC_mean"], cs_px.loc[BL, "MEC_lower"], cs_px.loc[BL, "MEC_upper"])
    else:
        mg = mr(a.loc[GSE, "MEC_mean"], a.loc[GSE, "MEC_lower"], a.loc[GSE, "MEC_upper"])
        mb = mr(a.loc[BL, "MEC_mean"], a.loc[BL, "MEC_lower"], a.loc[BL, "MEC_upper"])
    row = [lab, intc(a.loc[GSE, "median_nn_m"]), num(a.loc[GSE, "RMSE_mean"], 1), mg, mb]
    for j, reg in enumerate(REGIONS):
        v = bk[(bk.Design == d) & (bk.Scope == reg) & (bk.Model == GSE)].MEC_mean.iloc[0]
        row.append(num(v))
    rows.append(row)
hdr = [[{"text": "Design", "rowSpan": 2}, {"text": "Distance (m)", "rowSpan": 2},
        {"text": "All regions", "colSpan": 3}, {"text": "MEC of the GSE-only model in each region", "colSpan": 6}],
       ["RMSE, GSE-only (gC m⁻²)", "MEC, GSE-only", "MEC, baseline"] + REGIONS]
add("Table S8",
    "AGC accuracy of the GSE-only model and the mean-value baseline with spatial blocks of different "
    "size. The pixel hold-out is shown for comparison. Distance is the median distance from the test "
    "samples to the nearest training sample over all regions. Values are means of ten repeats, with the "
    "2.5th and 97.5th percentiles in parentheses.",
    hdr, rows, [12, 8, 10, 13, 13, 7, 7, 9.5, 7, 7, 7], landscape=True,
    csv_rows=[["Design", "Distance_m", "RMSE_GSE", "MEC_GSE", "MEC_baseline"] + REGIONS] + rows)

# ---------------------------------------------------------------------------
# Table S9  Each held-out survey and each held-out region
# ---------------------------------------------------------------------------
sv = f("agc_survey_transfer_ext.csv"); sv = sv[sv.Scope == "All other data"]
lo = f("agc_lolocv_pointTarget.csv")
dd = f("agc_design_distances_ext.csv")
def km(m):
    return num(m / 1000, 2) if m < 1000 else (num(m / 1000, 1) if m < 1e5 else intc(m / 1000))
rows = [[{"text": "(a) New survey, trained on all other data", "colSpan": 11}]]
for _, g in sv.groupby(["Region", "TestYear"], sort=True):
    t = f"{g.Region.iloc[0]} {g.TestYear.iloc[0]}"
    m = g.set_index("Model")
    dist = dd[(dd.Design == "Survey transfer") & (dd.Test == t)].median_nn_m.iloc[0]
    rows.append([t, intc(m.loc[GSE, "n_test"]), km(dist), num(m.loc[GSE, "mean_obs"], 1),
                 num(m.loc[PX, "MEC"]), num(m.loc[GSE, "MEC"]), num(m.loc[BL, "MEC"]),
                 num(m.loc[PX, "R2"]), num(m.loc[GSE, "R2"]),
                 pct(m.loc[PX, "Bias_pct"], 0), pct(m.loc[GSE, "Bias_pct"], 0), "–"])
rows.append([{"text": "(b) Unseen region, LOLO-CV", "colSpan": 11}])
for reg in REGIONS:
    m = lo[lo.TestLocation == reg].set_index("Model")
    dist = dd[(dd.Design == "Unseen regions (LOLO-CV)") & (dd.Test == reg)].median_nn_m.iloc[0]
    rows.append([reg, intc(m.loc[GSE, "n_test"]), km(dist), num(m.loc[GSE, "mean_obs"], 1),
                 num(m.loc[PX, "MEC"]), num(m.loc[GSE, "MEC"]), num(m.loc[BL, "MEC"]),
                 num(m.loc[PX, "R2"]), num(m.loc[GSE, "R2"]),
                 pct(m.loc[PX, "Bias_pct"], 0), pct(m.loc[GSE, "Bias_pct"], 0),
                 num(m.loc[PX, "R2_cindex_heldout"])])
# fix colSpan of the panel rows to the real number of columns
for r in rows:
    if isinstance(r[0], dict):
        r[0]["colSpan"] = 12
hdr = [[{"text": "Test survey or region", "rowSpan": 2}, {"text": "n", "rowSpan": 2},
        {"text": "Distance (km)", "rowSpan": 2}, {"text": "Mean AGC (gC m⁻²)", "rowSpan": 2},
        {"text": "MEC", "colSpan": 3}, {"text": "R²", "colSpan": 2}, {"text": "Bias (%)", "colSpan": 2},
        {"text": "Carbon index R²", "rowSpan": 2}],
       ["Proxy", "GSE-only", "Baseline", "Proxy", "GSE-only", "Proxy", "GSE-only"]]
csv9 = [["Test", "n", "Distance_km", "Mean_AGC", "MEC_proxy", "MEC_GSE", "MEC_baseline", "R2_proxy",
         "R2_GSE", "Bias_proxy", "Bias_GSE", "R2_cindex_heldout"]] + \
       [r for r in rows if not isinstance(r[0], dict)]
add("Table S9",
    "AGC accuracy for each held-out survey (new-survey design) and each held-out region (LOLO-CV). "
    "Distance is the median distance from the test samples to the nearest training sample. Mean AGC is "
    "the mean point estimate of the test samples. Bias is the bias of the mean AGC of the test samples. "
    "Carbon index R² is the accuracy of the carbon index model, retrained without the held-out "
    "region, in that region.",
    hdr, rows, [14, 7, 8, 9, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 9], landscape=True,
    note="Proxy, proxy workflow; baseline, mean-value baseline. R² is not defined for the baseline, "
         "which predicts one value for each test survey or region.", csv_rows=csv9)

# ---------------------------------------------------------------------------
# Table S10  Single-region transfer matrix (MEC)
# ---------------------------------------------------------------------------
cx = f("agc_cross_site_pointTarget.csv"); cx = cx[cx.Type == "Transfer"]
rows, bold = [], []
for mod, lab in [(PX, "(a) Proxy workflow"), (GSE, "(b) GSE-only model")]:
    rows.append([{"text": lab, "colSpan": 8}])
    x = cx[cx.Model == mod]
    for src in REGIONS:
        row = [src]
        for tst in REGIONS:
            if tst == src:
                row.append("–")
            else:
                v = x[(x.ModelSource == src) & (x.TestSite == tst)].MEC.iloc[0]
                row.append(num(v))
                if v > 0:
                    bold.append((len(rows), len(row) - 1))
        row.append(num(x[x.ModelSource == src].MEC.mean()))
        rows.append(row)
    rows.append(["Mean"] + [num(x[x.TestSite == t].MEC.mean()) for t in REGIONS] + [num(x.MEC.mean())])
    bold.append((len(rows) - 1, 0))
hdr = [[{"text": "Source region", "rowSpan": 2}, {"text": "Test region", "colSpan": 6},
        {"text": "Mean", "rowSpan": 2}], REGIONS]
add("Table S10",
    "MEC of single-region transfer. Each row is the source region in which the workflow was trained, "
    "and each column is the test region. Positive values are in bold. The last column and the last row "
    "of each panel give the means over the five test regions and the five source regions.",
    hdr, rows, [17, 11, 11, 14, 11, 11, 11, 12], bold=bold,
    csv_rows=[["Source"] + REGIONS + ["Mean"]] + [r for r in rows if not isinstance(r[0], dict)])

# ---------------------------------------------------------------------------
# Table S12  Uncertainty of the reference AGC values per region
# ---------------------------------------------------------------------------
ef = f("agc_error_floor_ext.csv").set_index("Region")
rows = [[r, intc(ef.loc[r, "n"]), num(ef.loc[r, "mean_AGC"], 1), num(ef.loc[r, "rmse_floor"], 1),
         num(100 * ef.loc[r, "rmse_floor"] / ef.loc[r, "mean_AGC"], 0),
         num(ef.loc[r, "RMSE_within_GSE"], 1)] for r in REGIONS + ["All regions"]]
add("Table S12",
    "Uncertainty of the reference AGC values in each region. For each sample, the 95% interval was "
    "converted to a standard deviation (SD) as (upper − lower)/3.92. The table gives the root mean "
    "square of these SDs. The RMSE of the GSE-only model in the pixel hold-out design is shown for "
    "comparison.",
    [["Region", "Samples (n)", "Mean AGC (gC m⁻²)", "Reference SD (gC m⁻²)",
      "SD as % of mean AGC", "RMSE, GSE-only (gC m⁻²)"]],
    rows, [20, 14, 16, 17, 16, 17], bold=[(6, c) for c in range(6)])

print("Tables written to", os.path.join(out_dir, "tables_csv"), ":", [t["id"] for t in tables])
