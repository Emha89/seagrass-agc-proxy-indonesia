"""
figure3_validation_designs.py
Figure 3: leak-free validation designs, drawn as Input / Process / Output.
Five designs ordered by distance (pixel, spatial blocks, new survey, unseen
region, single-region transfer), with the median distances of Table 3, the
cross-fitted proxy run and the sensitivity run with proxy models trained on all
samples. The quadrat sketch shows one transect with 3-5 quadrats per 10 m pixel.

Schematic only; no data are read. Region sample sizes are typed in from Table 2.
Usage (from the repository root):
  python3 R/figures_tables/figure3_validation_designs.py
Outputs (result/paper_tables/):
  Figure3_validation_design.png / .tiff (300 dpi), .pdf and .svg (editable text)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle

plt.rcParams["font.family"] = "Liberation Sans"   # metric-compatible with Arial

OUT_DIR = os.path.join("result", "paper_tables")
os.makedirs(OUT_DIR, exist_ok=True)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"             # keep text as editable text in SVG

W, H = 7.0, 6.45                                  # inches, full page width
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

# colours (Okabe-Ito, colour-blind safe)
TRAIN, TEST = "#0072B2", "#E69F00"
EDGE, FILL, BAND, SHARED = "#7F8794", "#F3F4F6", "#3B4A5A", "#E6EEF6"
TXT, GREY, LINE = "#1F2328", "#5A6270", "#B8BEC7"

regions = ["Ayau", "Bintan", "Karimunjawa", "Komodo", "Manado", "Rote"]
n_reg = [114, 59, 263, 1398, 292, 1550]


def box(x, y, w, h, fc=FILL, ec=EDGE, lw=0.7, r=0.04, z=1, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z, linestyle=ls,
                                boxstyle=f"round,pad=0,rounding_size={r}"))


def txt(x, y, s, size=7, weight="normal", color=TXT, ha="left", va="center",
        style="normal", z=5, **kw):
    ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha,
            va=va, fontstyle=style, zorder=z, **kw)


def arrow(x0, y0, x1, y1, color=GREY, lw=0.9, ms=7):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=4,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                mutation_scale=ms, shrinkA=0, shrinkB=0))


def cell(x, y, w, h, test):
    ax.add_patch(Rectangle((x, y), w, h, fc=TEST if test else TRAIN, ec="white", lw=0.6, zorder=3))


# ---------------------------------------------------------------- bands
X0, X1 = 0.42, 6.95
IN = (5.10, 6.40); PR = (1.00, 4.98); OUT = (0.05, 0.90)
for lab, (y0, y1) in (("INPUT", IN), ("PROCESS", PR), ("OUTPUT", OUT)):
    box(0.05, y0, 0.27, y1 - y0, fc=BAND, ec=BAND, r=0.03)
    txt(0.185, (y0 + y1) / 2, lab, size=8, weight="bold", color="white",
        ha="center", rotation=90)

# ---------------------------------------------------------------- INPUT
iy0, iy1 = IN
# I1: field samples and sample size per region
box(X0, iy0, 2.70, iy1 - iy0)
txt(X0 + 0.08, iy1 - 0.13, "Field samples", size=7.5, weight="bold")
txt(X0 + 0.08, iy1 - 0.29, "3,676 AGC quadrats in 1,767 pixels, one record per sample",
    size=6.5, color=GREY)
bx0, bmax = X0 + 0.80, 1.45
for i, (r, n) in enumerate(zip(regions, n_reg)):
    yy = iy1 - 0.47 - i * 0.135
    txt(bx0 - 0.05, yy, r, size=6.5, ha="right")
    ax.add_patch(Rectangle((bx0, yy - 0.045), bmax * n / max(n_reg), 0.09,
                           fc="#9AA3AF", ec="none", zorder=3))
    txt(bx0 + bmax * n / max(n_reg) + 0.04, yy, f"{n:,}", size=6.5)
txt(X0 + 2.62, iy1 - 0.47, "Proxy stages also use\nthe other samples\nin Table 2", size=6,
    color=GREY, ha="right", va="top", linespacing=1.25)

# I2: photo quadrats along one diagonal transect inside 10-m pixels
x2 = X0 + 2.83
box(x2, iy0, 1.62, iy1 - iy0)
txt(x2 + 0.08, iy1 - 0.13, "Quadrats and pixels", size=7.5, weight="bold")
cs, gx, gy = 0.20, x2 + 0.10, iy0 + 0.27
for c in range(3):
    for r in range(3):
        ax.add_patch(Rectangle((gx + c * cs, gy + r * cs), cs, cs, fc="white", ec=EDGE, lw=0.6, zorder=3))
# one diagonal transect; 4 photo quadrats in each 10-m pixel it crosses
ax.plot([gx + 0.08 * cs, gx + 2.92 * cs], [gy + 0.08 * cs, gy + 2.92 * cs], color=LINE, lw=0.8, zorder=3.5)
for c in range(3):
    for f in (0.2, 0.4, 0.6, 0.8):
        ax.add_patch(Circle((gx + (c + f) * cs, gy + (c + f) * cs), 0.015, fc=TXT, ec="none", zorder=4))
txt(gx + 1.5 * cs, gy - 0.10, "10-m pixels", size=6, color=GREY, ha="center")
txt(gx + 3 * cs + 0.08, gy + 1.5 * cs,
    "3–5 photo quadrats\nper pixel along a\ntransect; 71% of\nquadrats share a\npixel and identical\nGSE values",
    size=6.1, linespacing=1.25)

# I3: predictors and target
x3 = x2 + 1.75
box(x3, iy0, X1 - x3, iy1 - iy0)
txt(x3 + 0.08, iy1 - 0.13, "Predictors", size=7.5, weight="bold")
txt(x3 + 0.08, iy1 - 0.42, "64 GSE bands (annual, 10 m)\nWater depth\nDistance to land",
    size=6.5, linespacing=1.35)
txt(x3 + 0.08, iy0 + 0.40, "Target", size=7.5, weight="bold")
txt(x3 + 0.08, iy0 + 0.19, "AGC point estimate (gC m$^{-2}$)", size=6.5)

# ---------------------------------------------------------------- PROCESS
py1 = PR[1]
txt(X0, py1 - 0.13, "Split of field samples", size=7.5, weight="bold")
lx = X0 + 1.55
for lab, col in (("Training", TRAIN), ("Test", TEST)):
    ax.add_patch(Rectangle((lx, py1 - 0.18), 0.16, 0.10, fc=col, ec="none", zorder=3))
    txt(lx + 0.21, py1 - 0.13, lab, size=6.5)
    lx += 0.80
txt(X0 + 4.55, py1 - 0.13, "ordered by distance to training data", size=6.3, color=GREY,
    style="italic", ha="right")

barx0, barx1 = X0 + 0.07, X0 + 3.80
bw = (barx1 - barx0) / 6
for i, r in enumerate(regions):
    txt(barx0 + (i + 0.5) * bw, py1 - 0.34, r, size=6.3, ha="center", color=GREY)

lanes = [
    ("A", "Pixel hold-out", "median distance 6 m", "×10 repeats",
     "30% of pixels in each region held out; quadrats of one pixel kept together"),
    ("B", "Spatial blocks", "median distance 43–663 m", "×10 repeats",
     "30% of 100, 250, 500 or 1000 m blocks in each region held out"),
    ("C", "New survey", "median distance 57 m–23 km", "×8 surveys",
     "One survey year of one region held out (one shown); trained on all other data"),
    ("D", "Unseen region (LOLO-CV)", "median distance 411–976 km", "×6 folds",
     "Each region held out once (one fold shown)"),
    ("E", "Single-region transfer", "median distance 411–2,933 km", "×30 pairs",
     "Each region used once as source; each other region scored separately"),
]
pix_test = {0: [2, 6, 8], 1: [1, 4, 7], 2: [0, 5, 8], 3: [3, 6, 9], 4: [1, 2, 7], 5: [4, 5, 8]}
blk_test = {0: 2, 1: 0, 2: 3, 3: 1, 4: 2, 5: 0}
lane_w, lane_h, gap_l = 4.55, 0.58, 0.05
lane_top = py1 - 0.43
lane_y = []
for k, (tag, title, dist, rep, note) in enumerate(lanes):
    top = lane_top - k * (lane_h + gap_l)
    box(X0, top - lane_h, lane_w, lane_h, fc="white", ec=EDGE, lw=0.6)
    txt(X0 + 0.07, top - 0.11, tag, size=8, weight="bold")
    txt(X0 + 0.23, top - 0.11, title, size=7, weight="bold")
    txt(X0 + 2.05, top - 0.11, dist, size=6.5, color=GREY, style="italic")
    by, bh = top - 0.35, 0.15
    for i in range(6):
        x = barx0 + i * bw
        if tag == "A":                              # 10 pixels per region, 3 held out
            pw = bw / 10
            for j in range(10):
                cell(x + j * pw, by, pw, bh, j in pix_test[i])
        elif tag == "B":                            # 4 blocks per region, 1 held out
            qw = bw / 4
            for j in range(4):
                cell(x + j * qw, by, qw, bh, j == blk_test[i])
        elif tag == "C":                            # surveys of Karimunjawa, one held out
            if i == 2:
                sw = bw / 3
                for j in range(3):
                    cell(x + j * sw, by, sw, bh, j == 0)
                txt(x + sw / 2, by + bh / 2, "2018", size=5, color=TXT, ha="center")
            else:
                cell(x, by, bw, bh, False)
        else:
            is_test = (i == 4) if tag == "D" else (i != 3)
            cell(x, by, bw, bh, is_test)
    txt(barx1 + 0.08, by + bh / 2, rep, size=7, weight="bold")
    txt(barx0, by - 0.10, note, size=6, color=GREY)
    lane_y.append(by + bh / 2)
lanes_bot = lane_top - 5 * lane_h - 4 * gap_l

# shared box: steps in every repeat, fold, survey and pair
sx0, sx1 = X0 + lane_w + 0.15, X1
s_top, s_bot = py1 - 0.02, lanes_bot
box(sx0, s_bot, sx1 - sx0, s_top - s_bot, fc=SHARED, ec=EDGE)
cx = (sx0 + sx1) / 2
txt(cx, s_top - 0.15, "In every repeat, fold,\nsurvey and pair", size=7, weight="bold",
    ha="center", linespacing=1.15)
steps = [("Training samples only", None, "-"),
         ("Refit the 5 proxy models", "presence, morphology, cover,\nAGB, carbon index", "-"),
         ("Cross-fitted proxy values", "extra run in A and B only:\nout-of-fold values for training samples", "--"),
         ("Fit the AGC models", None, "-"),
         ("Predict the test samples", None, "-")]
sy = s_top - 0.40
for j, (s, sub, ls) in enumerate(steps):
    hh = 0.34 if sub else 0.19
    box(sx0 + 0.10, sy - hh, sx1 - sx0 - 0.20, hh, fc="white", ec=EDGE, lw=0.6, r=0.03, z=2, ls=ls)
    if sub:
        txt(cx, sy - 0.10, s, size=6.5, weight="bold", ha="center")
        txt(cx, sy - 0.24, sub, size=6, color=GREY, ha="center", linespacing=1.15)
    else:
        txt(cx, sy - hh / 2, s, size=6.5, weight="bold", ha="center")
    if j < len(steps) - 1:
        arrow(cx, sy - hh, cx, sy - hh - 0.08, ms=5)
    sy -= hh + 0.08
txt(sx0 + 0.12, sy - 0.02, "Models compared on\nthe same test samples",
    size=6.3, weight="bold", va="top", linespacing=1.2)
txt(sx0 + 0.12, sy - 0.32,
    "• Proxy workflow\n• Proxy workflow, cross-fitted\n   (A and B)\n• GSE-only model\n• Mean-value baseline\n   (training mean)",
    size=6.2, va="top", linespacing=1.3)
for ym in lane_y:                                  # lanes -> shared box
    arrow(X0 + lane_w, ym, sx0, ym, ms=6)

# sensitivity lane
d_top, d_bot = lanes_bot - 0.08, PR[0] + 0.05
box(X0, d_bot, X1 - X0, d_top - d_bot, fc="white", ec=EDGE, lw=0.6, ls="--")
txt(X0 + 0.07, (d_top + d_bot) / 2, "Sensitivity", size=7, weight="bold")
txt(X0 + 0.85, (d_top + d_bot) / 2,
    "Design A repeated with the proxy models trained on all samples, including the test samples,"
    " to show how much refitting inside each split matters", size=6.3, color=GREY)

# ---------------------------------------------------------------- OUTPUT
outs = [("Metrics", "RMSE, MAE, R², MEC and\nregional mean bias (%),\nscored against the field\npoint estimates"),
        ("A and B", "Mean and 2.5th–97.5th\npercentiles over 10 repeats,\noverall and per region;\npaired test of RMSE"),
        ("C, D and E", "Mean and range over\n8 surveys, 6 regions\nor 30 region pairs,\nequal weight"),
        ("Distance", "Median distance to the\nnearest training sample;\nmapped area by distance\nto the nearest AGC sample"),
        ("Sensitivity", "Apparent accuracy with\nproxy models trained\non all samples")]
oy0, oy1 = OUT
og = 0.125
ow = (X1 - X0 - 4 * og) / 5
for k, (head, body) in enumerate(outs):
    x = X0 + k * (ow + og)
    box(x, oy0, ow, oy1 - oy0, fc=SHARED if k == 0 else FILL)
    txt(x + 0.08, oy1 - 0.12, head, size=7.5, weight="bold")
    txt(x + 0.08, oy1 - 0.25, body, size=6.2, va="top", linespacing=1.3)
    arrow(x + ow / 2, PR[0], x + ow / 2, oy1)

# input -> process arrows
arrow(X0 + 1.35, iy0, X0 + 1.35, py1)
arrow(x2 + 0.81, iy0, x2 + 0.81, py1)
arrow(x3 + (X1 - x3) / 2, iy0, x3 + (X1 - x3) / 2, s_top)

for ext, kw in (("png", {"dpi": 300}), ("tiff", {"dpi": 300, "pil_kwargs": {"compression": "tiff_lzw"}}),
                ("pdf", {}), ("svg", {})):
    fig.savefig(os.path.join(OUT_DIR, f"Figure3_validation_design.{ext}"), **kw)
print("done")
