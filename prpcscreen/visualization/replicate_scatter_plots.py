from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def _candidate_mask(df: pd.DataFrame) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    nt = df.get("Is_NT_ctrl", pd.Series(False, index=df.index)).astype(bool)
    pos = df.get("Is_pos_ctrl", pd.Series(False, index=df.index)).astype(bool)
    is_control = nt | pos

    if "Gene_symbol" in df.columns:
        gene = df["Gene_symbol"].astype(str).str.strip()
        gene_ok = df["Gene_symbol"].notna() & gene.ne("")
        if "Entrez_ID" in df.columns:
            entrez_ok = pd.to_numeric(df["Entrez_ID"], errors="coerce").notna()
            mask &= is_control | gene_ok | entrez_ok
        else:
            mask &= is_control | gene_ok
    elif "Entrez_ID" in df.columns:
        mask &= is_control | pd.to_numeric(df["Entrez_ID"], errors="coerce").notna()

    return mask


def replicate_diagnostics(df: pd.DataFrame, stem: str = "Raw"):
    keep = _candidate_mask(df)
    x = pd.to_numeric(df[f"{stem}_rep1"], errors="coerce")
    y = pd.to_numeric(df[f"{stem}_rep2"], errors="coerce")
    valid = keep & x.notna() & y.notna()
    removed = int((~keep).sum())

    if valid.sum() == 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f"No valid {stem} replicate values found", ha="center", va="center")
        ax.set_axis_off()
        return fig, ax

    is_nt = df["Is_NT_ctrl"].astype(bool) if "Is_NT_ctrl" in df else pd.Series(False, index=df.index)
    is_pos = df["Is_pos_ctrl"].astype(bool) if "Is_pos_ctrl" in df else pd.Series(False, index=df.index)
    is_gene = (~is_nt) & (~is_pos)

    xv = x[valid]
    yv = y[valid]
    mv = (xv + yv) / 2.0
    dv = yv - xv
    av = mv.abs()
    adv = dv.abs()

    gene_mask = (valid & is_gene)[valid]
    nt_mask = (valid & is_nt)[valid]
    pos_mask = (valid & is_pos)[valid]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    ax1, ax2, ax3, ax4 = axes.flatten()

    # Panel 1: replicate scatter with controls highlighted.
    ax1.scatter(xv[gene_mask], yv[gene_mask], s=8, alpha=0.25, c="#7a7a7a", label=f"Genes (n={int(gene_mask.sum())})")
    ax1.scatter(xv[nt_mask], yv[nt_mask], s=12, alpha=0.9, c="#d62728", label=f"Non-targeting (n={int(nt_mask.sum())})")
    ax1.scatter(xv[pos_mask], yv[pos_mask], s=12, alpha=0.9, c="#1f77b4", label=f"Positive (n={int(pos_mask.sum())})")
    lo = float(np.nanpercentile(pd.concat([xv, yv]), 0.5))
    hi = float(np.nanpercentile(pd.concat([xv, yv]), 99.5))
    if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
        pad = (hi - lo) * 0.05
        ax1.set_xlim(lo - pad, hi + pad)
        ax1.set_ylim(lo - pad, hi + pad)
    ax1.plot(ax1.get_xlim(), ax1.get_xlim(), ls="--", lw=1, c="black", alpha=0.7)
    r = np.corrcoef(xv, yv)[0, 1] if xv.size > 1 else np.nan
    ax1.text(0.02, 0.98, f"r={r:.3f}", transform=ax1.transAxes, va="top")
    ax1.set_title("Replicate Scatter")
    ax1.set_xlabel(f"{stem}_rep1")
    ax1.set_ylabel(f"{stem}_rep2")
    if removed:
        ax1.text(0.99, 0.01, f"Filtered rows: {removed}", transform=ax1.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")
    ax1.legend(loc="lower right", frameon=False, fontsize=8)

    # Panel 2: Bland-Altman (difference vs mean).
    ax2.scatter(mv[gene_mask], dv[gene_mask], s=8, alpha=0.25, c="#7a7a7a")
    ax2.scatter(mv[nt_mask], dv[nt_mask], s=12, alpha=0.9, c="#d62728")
    ax2.scatter(mv[pos_mask], dv[pos_mask], s=12, alpha=0.9, c="#1f77b4")
    bias = float(np.nanmean(dv))
    sd = float(np.nanstd(dv))
    ax2.axhline(bias, color="black", lw=1)
    if np.isfinite(sd):
        ax2.axhline(bias + 1.96 * sd, color="black", lw=1, ls=":")
        ax2.axhline(bias - 1.96 * sd, color="black", lw=1, ls=":")
    ax2.set_title("Bland-Altman")
    ax2.set_xlabel("Mean effect")
    ax2.set_ylabel("Rep2 - Rep1")

    # Panel 3: absolute disagreement as a function of absolute effect size.
    ax3.scatter(av[gene_mask], adv[gene_mask], s=8, alpha=0.2, c="#7a7a7a")
    ax3.scatter(av[nt_mask], adv[nt_mask], s=12, alpha=0.8, c="#d62728")
    ax3.scatter(av[pos_mask], adv[pos_mask], s=12, alpha=0.8, c="#1f77b4")
    if av.size >= 100:
        bins = pd.qcut(av, q=min(20, max(5, int(np.sqrt(av.size)))), duplicates="drop")
        trend = pd.DataFrame({"abs_mean": av, "abs_delta": adv, "bin": bins}).groupby("bin", observed=True).median(numeric_only=True)
        ax3.plot(trend["abs_mean"], trend["abs_delta"], color="black", lw=1.5)
    ax3.set_title("Error vs Effect Magnitude")
    ax3.set_xlabel("|Mean effect|")
    ax3.set_ylabel("|Rep2 - Rep1|")

    # Panel 4: binned replicate correlation vs |mean effect|.
    if av.size >= 60:
        bin_codes = pd.qcut(av, q=min(10, max(5, int(np.sqrt(av.size) / 2))), duplicates="drop")
        tmp = pd.DataFrame({"x": xv, "y": yv, "abs_mean": av, "bin": bin_codes})
        rows: list[tuple[float, float, int]] = []
        for _, sub in tmp.groupby("bin", observed=True):
            n = len(sub)
            if n < 30:
                continue
            sx = float(sub["x"].std())
            sy = float(sub["y"].std())
            if np.isclose(sx, 0) or np.isclose(sy, 0):
                continue
            corr = float(sub["x"].corr(sub["y"]))
            center = float(sub["abs_mean"].median())
            rows.append((center, corr, n))
        if rows:
            corr_df = pd.DataFrame(rows, columns=["center", "corr", "n"]).sort_values("center")
            ax4.plot(corr_df["center"], corr_df["corr"], marker="o", lw=1.5, c="#2f4b7c")
            for _, row in corr_df.iterrows():
                ax4.text(row["center"], row["corr"], f"n={int(row['n'])}", fontsize=7, ha="left", va="bottom")
    ax4.axhline(0, lw=1, ls=":", color="grey")
    ax4.set_ylim(-1.0, 1.0)
    ax4.set_title("Binned Correlation vs |Effect|")
    ax4.set_xlabel("|Mean effect|")
    ax4.set_ylabel("corr(rep1, rep2)")

    legend_text = (
        "Interpretation legend\n"
        "Panel 1 (Replicate Scatter): each point is one perturbation with x=rep1 and y=rep2.\n"
        "Good: points close to identity line y=x, Pearson r close to 1, controls behaving as expected.\n"
        "Rule of thumb: r >= 0.70 acceptable, >= 0.80 good, >= 0.90 excellent (screen/platform dependent).\n"
        "Panel 2 (Bland-Altman): x=(rep1+rep2)/2, y=(rep2-rep1). Center line is bias=mean(delta).\n"
        "Dotted lines are limits of agreement: bias +/- 1.96*SD(delta).\n"
        "Good: bias near 0 and most points within LoA without strong trends vs mean effect.\n"
        "Panel 3 (Error vs Effect): x=|mean|, y=|delta| with median trend line.\n"
        "Good: flat or gently rising trend; steep rise indicates heteroscedastic noise at strong effects.\n"
        "Panel 4 (Binned Correlation): corr(rep1,rep2) within bins of |mean effect|.\n"
        "Good: correlation remains positive and fairly stable across bins; drops toward 0/negative are concerning.\n"
        "Quality context: thresholds are guidance only. Judge alongside control separation, hit reproducibility, and assay biology."
    )

    fig.suptitle(f"Replicate Diagnostics ({stem})", y=0.98)
    fig.text(
        0.01,
        0.01,
        legend_text,
        ha="left",
        va="bottom",
        fontsize=8,
        color="#333333",
        linespacing=1.15,
    )
    fig.tight_layout(rect=(0.0, 0.22, 1.0, 0.96))
    return fig, axes
