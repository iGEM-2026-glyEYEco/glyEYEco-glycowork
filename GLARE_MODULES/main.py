"""

DOCUMENTATION FOR GLARE (Glycan-Lectin Analysis for Recognition and Evaluation)

ALL COMBINATIONS FOR EXAMPLES REGARDING TERMINAL USAGE
python main.py --help
python main.py --list-systems
python main.py --cell eye
python main.py --mucin MUC16
python main.py --mucin-map
python main.py --compare eye gut_lining
python main.py --motif core1 eye gut_lining
python main.py --binding Gal-3
python main.py --binding-pair "Gal-3" "core1"
python main.py --binding-galectins eye
python main.py --tissue-rank eye
python main.py --binding-rank core1
python main.py --off-target eye
python main.py --motif-affinity eye
python main.py --target-chart eye
python main.py --disease-compare dry_eye
python main.py --disease-compare vaginal_cancer
python main.py --disease-compare kidney_cancer
python main.py --disease-compare lung_cancer
python main.py --cross-system-synthesis

"""

from __future__ import annotations

import argparse
import functools
import re
import sys
import textwrap
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    from rich.console import Console
    from rich.table   import Table
    from rich.panel   import Panel
    from rich.rule    import Rule
    from rich.prompt  import Prompt
    from rich         import box as rbox
    HAS_RICH = True
except ImportError:
    HAS_RICH = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

try:
    from scipy import stats as _scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from GLARE_MODULES.databases.glycan_data import (GLYCAN_TYPES, SYSTEMS, TISSUES, MUCINS, TISSUE_TO_MUCINS, MOTIFS, GLYCAN_ALIASES, TARGET_CLASSES, DISEASE_SIGNATURES)
from GLARE_MODULES.databases.lectin_data import LECTIN_DATABASE

CITATIONS_BASE = [
    ("glycowork",  "Thomès et al. 2021 — glycowork: a Python package for glycoinformatics. "
                   "Cell Rep. https://doi.org/10.1016/j.celrep.2021.109819"),
    ("GlyTouCan", "Fujita et al. 2021 — GlyTouCan: The international glycan structure repository. "
                  "https://doi.org/10.1093/nar/gkab918"),
    ("GlyCosMos", "https://glycosmos.org/"),
]

CITATIONS_TISSUE: dict[str, list] = {
    "eye": [
        ("Argüeso 2008",  "Mantelli & Argüeso — Glycoconjugate composition of the ocular surface. "
                          "https://doi.org/10.1016/j.preteyeres.2008.05.002"),
        ("Royle 2003",    "Royle et al. — Sialylated N-glycans on ocular glycoproteins. "
                          "https://doi.org/10.1074/jbc.M213003200"),
    ],
}

MOTIF_CITATIONS = [
    ("Taylor 2015",  "Taylor & Drickamer — Paradigms for glycan-binding receptors in cell adhesion. "
                     "https://doi.org/10.1016/j.cub.2015.03.022"),
    ("Varki 2017",   "Varki et al. — Biological roles of glycans. "
                     "https://doi.org/10.1091/mbc.e16-10-0739"),
]

LECTIN_CITATIONS = [
    ("Schauer 2009",   "Schauer — Sialic acids as regulators of molecular and cellular interactions. "
                       "https://doi.org/10.1016/j.cub.2009.08.016"),
    ("Drickamer 2015", "Drickamer & Taylor — Recent insights into structures and functions of "
                       "C-type lectins. https://doi.org/10.1016/j.cub.2015.09.051"),
]

MUCIN_CITATIONS = [
    ("Bansil 2006",        "Bansil & Turner — Mucin structure, aggregation, physiological functions. "
                           "https://doi.org/10.1016/j.cocis.2006.04.001"),
    ("Hollingsworth 2004", "Hollingsworth & Swanson — Mucins in cancer. "
                           "https://doi.org/10.1038/nrc1251"),
]

GALECTIN_CITATIONS = [
    ("Johannes 2018", "Johannes et al. — Galectins at a glance. "
                      "https://doi.org/10.1242/jcs.208884"),
    ("Leffler 2002",  "Leffler et al. — The human galectin gene family. "
                      "https://doi.org/10.1074/jbc.R200020200"),
]

@functools.lru_cache(maxsize=1)
def load_data():
    try:
        from glycowork.glycan_data.loader import df_glycan, glycan_binding
        return df_glycan, glycan_binding
    except ImportError:
        print("[ERROR] glycowork not installed.  Run:  pip install glycowork")
        sys.exit(1)

@functools.lru_cache(maxsize=1)
def load_lectin_lib():
    from glycowork.motif.analysis import load_lectin_lib as _load
    return _load()

try:
    from rich.console import Console
    console = Console(highlight=False)
except ImportError:
    console = None

OUTPUT_DIR = Path("output")
CAT_COLORS = ["#4878CF", "#D65F00", "#3A9E4A", "#B5282F",
               "#8862AE", "#987845", "#D16DA8"]
GTYPE_VIZ = {"N": CAT_COLORS[0], "O": CAT_COLORS[1],
              "free": CAT_COLORS[2], "lipid": CAT_COLORS[3], "repeat": CAT_COLORS[4]}


def _ensure_output():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR


def rprint(msg: str = ""):
    if HAS_RICH:
        console.print(msg)
    else:
        print(re.sub(r"\[/?[^\]]*\]", "", msg))


def rule(title: str = ""):
    if HAS_RICH:
        console.print(Rule(f" {title} ", style="bold cyan"))
    else:
        print("\n" + "─" * 2 + f"  {title}  " + "─" * 60)


def panel(body: str, title: str = "", style: str = "cyan"):
    if HAS_RICH:
        console.print(Panel(body, title=title, border_style=style, padding=(0, 1)))
    else:
        if title:
            print(f"\n── {title} ──")
        print(re.sub(r"\[/?[^\]]*\]", "", body))


# TISSUE FUNCTION HELPERS
def col_contains_any(val: Any, terms: list[str]) -> bool:
    if val is None or (isinstance(val, float)):
        return False
    items = val if isinstance(val, list) else [str(val)]
    return any(t.lower() in item.lower() for item in items for t in terms)


def get_tissue_glycans(tissue_key: str):
    """
    Universal tissue filter. Returns df with added 'organ' column.
    """
    df_glycan, _ = load_data()
    info  = TISSUES[tissue_key]
    terms = info["terms"]
    mask  = df_glycan["tissue_sample"].apply(lambda v: col_contains_any(v, terms))
    df    = df_glycan[mask].copy().reset_index(drop=True)
    df["organ"] = info["display"]
    return df


def tissue_keys_from_query(query: str) -> list[str]:
    q = query.lower().replace("-", "_").replace(" ", "_")
    if q in TISSUES:
        return [q]
    for k, info in TISSUES.items():
        if q == info["display"].lower().replace(" ", "_"):
            return [k]
    matches = [k for k in TISSUES if q in k]
    if matches:
        return matches
    return []


# MUCIN FUNCTION HELPERS
def resolve_mucin(query: str) -> tuple[str | None, dict | None]:
    def norm(s): return re.sub(r"[\s\-_]", "", s.lower())
    q, qn = query.lower().strip(), norm(query)
    for key, info in MUCINS.items():
        if q in info["aliases"] or qn == norm(key) or any(qn == norm(a) for a in info["aliases"]):
            return key, info
    for key, info in MUCINS.items():
        if qn in norm(key) or any(qn in norm(a) for a in info["aliases"]):
            return key, info
    return None, None


def resolve_glycan(query: str) -> tuple[str, str]:
    q = query.strip().lower()
    if q in GLYCAN_ALIASES:
        iupac = GLYCAN_ALIASES[q]
        return iupac, f"{query} -> {iupac}"
    q2 = re.sub(r"\s+", " ", q)
    if q2 in GLYCAN_ALIASES:
        iupac = GLYCAN_ALIASES[q2]
        return iupac, f"{query} -> {iupac}"
    return query.strip(), query.strip()


@functools.lru_cache(maxsize=1)
def _build_human_lectin_index() -> dict[str, dict]:
    idx = {}
    for key, info in LECTIN_DATABASE.items():
        idx[key.lower()] = info
        idx[re.sub(r"[\s\-_]", "", key.lower())] = info
        for a in info["abbreviation"]:
            idx[a.lower()] = info
            idx[re.sub(r"[\s\-_]", "", a.lower())] = info
        idx[info["gene"].lower()] = info
    return idx


def resolve_lectin(query: str) -> dict | None:
    idx = _build_human_lectin_index()
    q  = query.lower().strip()
    qn = re.sub(r"[\s\-_]", "", q)
    if q  in idx: return idx[q]
    if qn in idx: return idx[qn]
    for k, v in idx.items():
        if qn in k or k in qn:
            return v
    return None


def get_all_galectins() -> list[dict]:
    return [info for key, info in LECTIN_DATABASE.items()
            if info["family"] == "Galectin"]


def _fig_style(ax, title="", xlabel="", ylabel=""):
    ax.set_title(title, fontsize=11, pad=8, fontweight="bold", color="#1a1a1a")
    if xlabel: ax.set_xlabel(xlabel, fontsize=9, color="#555")
    if ylabel: ax.set_ylabel(ylabel, fontsize=9, color="#555")
    for spine in ("top", "right"): ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"): ax.spines[spine].set_color("#ddd")
    ax.tick_params(colors="#555", labelsize=8)
    ax.yaxis.grid(True, color="#eee", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def _save(fig, name: str) -> Path:
    d  = _ensure_output()
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    p  = d / f"{name}_{ts}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    rprint(f"  [dim]Chart -> {p}[/dim]")
    return p


def save_csv(df, name: str) -> Path:
    d  = _ensure_output()
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    p  = d / f"{name}_{ts}.csv"
    df.to_csv(p, index=False)
    rprint(f"  [dim]CSV  -> {p}[/dim]")
    return p


def plot_gtype_comp_bar(df, title: str, tag: str) -> Path | None:
    if not HAS_PLOT: return None
    counts = df["glycan_type"].value_counts(dropna=True)
    total  = counts.sum()
    if total == 0: return None
    labels = [GLYCAN_TYPES.get(k, k) for k in counts.index]
    fracs  = counts.values / total
    colors = [GTYPE_VIZ.get(k, "#888") for k in counts.index]
    fig, ax = plt.subplots(figsize=(8, 1.8))
    left = 0.0
    for frac, color, label in zip(fracs, colors, labels):
        ax.barh(0, frac, left=left, color=color, height=0.5, edgecolor="white", linewidth=0.8)
        if frac > 0.04:
            ax.text(left + frac / 2, 0, f"{frac*100:.0f}%",
                    va="center", ha="center", fontsize=8, color="white", fontweight="bold")
        left += frac
    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, labels)]
    ax.legend(handles=patches, loc="upper center", bbox_to_anchor=(0.5, -0.35),
              ncol=min(len(labels), 4), fontsize=8, frameon=False)
    ax.set_xlim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=8)
    ax.set_title(title, fontsize=11, pad=8, fontweight="bold")
    for spine in ("top", "left", "right"): ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#ddd")
    plt.tight_layout()
    return _save(fig, tag + "_comp")


def plot_gtype_bar(df, title: str, tag: str) -> Path | None:
    if not HAS_PLOT: return None
    counts = df["glycan_type"].value_counts(dropna=True)
    n      = len(df)
    labels = [GLYCAN_TYPES.get(k, k) for k in counts.index]
    colors = [GTYPE_VIZ.get(k, "#888") for k in counts.index]
    fig, ax = plt.subplots(figsize=(7, max(3, len(labels) * 0.55 + 1.2)))
    bars = ax.barh(labels, counts.values, color=colors, height=0.6, edgecolor="white")
    for bar, val in zip(bars, counts.values):
        ax.text(val + counts.values.max() * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val}  ({100*val/n:.0f}%)", va="center", fontsize=8, color="#333")
    ax.invert_yaxis()
    ax.set_xlim(0, counts.values.max() * 1.28)
    _fig_style(ax, title=title, xlabel="Number of glycans")
    plt.tight_layout()
    return _save(fig, tag + "_bar")


def plot_motif_bar(motif_hits: dict, title: str, tag: str, total: int = 0) -> Path | None:
    if not HAS_PLOT or not motif_hits: return None
    items  = sorted(motif_hits.items(), key=lambda x: -x[1])
    labels = [i[0] for i in items]
    values = [i[1] for i in items]
    colors = [CAT_COLORS[i % len(CAT_COLORS)] for i in range(len(labels))]
    fig, ax = plt.subplots(figsize=(8, max(3, len(labels) * 0.5 + 1.5)))
    bars = ax.barh(labels, values, color=colors, height=0.6, edgecolor="white")
    for bar, val in zip(bars, values):
        pct = f"  {val}" if not total else f"  {val} ({100*val/total:.0f}%)"
        ax.text(val + max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
                pct, va="center", fontsize=8, color="#333")
    ax.invert_yaxis()
    ax.set_xlim(0, max(values) * 1.3)
    _fig_style(ax, title=title, xlabel="Glycans containing motif")
    plt.tight_layout()
    return _save(fig, tag + "_motifs")


def plot_compare_grouped(data: dict, title: str, tag: str) -> Path | None:
    if not HAS_PLOT: return None
    all_types = sorted({t for s in data.values() for t in s.index})
    x     = np.arange(len(all_types))
    n     = len(data)
    width = 0.7 / n
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (label, series) in enumerate(data.items()):
        vals = [series.get(t, 0) for t in all_types]
        ax.bar(x + (i - (n-1)/2) * width, vals, width * 0.9,
               label=label, color=CAT_COLORS[i % len(CAT_COLORS)], edgecolor="white")
    ax.set_xticks(x)
    ax.set_xticklabels([GLYCAN_TYPES.get(t, t) for t in all_types], rotation=25, ha="right")
    ax.legend(fontsize=8, frameon=False)
    _fig_style(ax, title=title, ylabel="Count")
    plt.tight_layout()
    return _save(fig, tag)


def plot_binding_hbar(scores: list[tuple], title: str, tag: str) -> Path | None:
    if not HAS_PLOT or not scores: return None
    labels = [s[0] for s in scores]
    values = [s[1] for s in scores]
    cats   = [s[2] for s in scores]
    cmap   = {"Primary": CAT_COLORS[2], "Secondary": CAT_COLORS[1], "None": CAT_COLORS[3]}
    colors = [cmap.get(c, "#888") for c in cats]
    fig, ax = plt.subplots(figsize=(8, max(3, len(labels) * 0.5 + 1.5)))
    ax.barh(labels, values, color=colors, height=0.6, edgecolor="white")
    for bar, val, cat in zip(ax.patches, values, cats):
        ax.text(max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
                f"  {cat}  ({val:.2f})", va="center", fontsize=8, color="#222")
    patches = [mpatches.Patch(color=v, label=k) for k, v in cmap.items() if k != "None"]
    ax.legend(handles=patches, fontsize=8, frameon=False, loc="lower right")
    ax.invert_yaxis()
    _fig_style(ax, title=title, xlabel="Binding level (0=none, 1=primary, 2=secondary)")
    plt.tight_layout()
    return _save(fig, tag)


def print_citations(extras: list | None = None, tissue_key: str | None = None):
    rprint("\n[bold yellow]Citations & data sources[/bold yellow]")
    all_cites = list(CITATIONS_BASE)
    if tissue_key and tissue_key in CITATIONS_TISSUE:
        all_cites = CITATIONS_TISSUE[tissue_key] + all_cites
    if extras:
        all_cites = extras + all_cites
    for ref, desc in all_cites:
        rprint(f"  [cyan]{ref}[/cyan]  {desc}")


def cmd_list_systems():
    rule("Biological Systems")
    for sys_name, tkeys in SYSTEMS:
        rprint(f"\n[bold yellow]{sys_name}[/bold yellow]")
        for tk in tkeys:
            info   = TISSUES[tk]
            mucins = ", ".join(TISSUE_TO_MUCINS.get(tk, [])) or "—"
            rprint(f"  [cyan]• {info['display']}[/cyan]  [dim](key: {tk})[/dim]")
            rprint(f"    [dim]{info['description'][:120]}...[/dim]")
            if TISSUE_TO_MUCINS.get(tk):
                rprint(f"    [dim]Mucins: {mucins}[/dim]")


def cmd_list_tissues():
    rule("Tissues")
    for sys_name, tkeys in SYSTEMS:
        rprint(f"\n[bold yellow]{sys_name}[/bold yellow]")
        for tk in tkeys:
            info   = TISSUES[tk]
            mucins = ", ".join(TISSUE_TO_MUCINS.get(tk, [])) or "—"
            is_sub = info.get("parent") is not None
            indent = "    " if is_sub else "  "
            prefix = "╰─ " if is_sub else "◆ "
            star   = " [bold red]★[/bold red]" if info.get("cell_type") else ""
            if HAS_RICH:
                rprint(
                    f"{indent}[cyan]{prefix}{tk}[/cyan]"
                    f"  [white]{info['display']}[/white]{star}"
                    f"  [dim]{mucins}[/dim]"
                )
                if info.get("primary_glycan"):
                    rprint(f"      [green dim]Primary glycan: {info['primary_glycan']}[/green dim]")
                if info.get("primary_lectin"):
                    rprint(f"      [green dim]Primary lectin: {info['primary_lectin']}[/green dim]")
            else:
                print(f"{indent}{tk:28s}  {info['display']}")


def cmd_list_mucins():
    rule("Mucin Knowledge Base")
    if HAS_RICH:
        tbl = Table(box=rbox.SIMPLE_HEAVY)
        tbl.add_column("Name",         style="cyan",   no_wrap=True)
        tbl.add_column("Type",         style="yellow")
        tbl.add_column("Key tissues")
        tbl.add_column("Glycan types", style="green")
        for name, info in MUCINS.items():
            ts = "; ".join(info["tissues"][:2])
            if len(info["tissues"]) > 2: ts += " …"
            tbl.add_row(name, info["type"], ts, ", ".join(info["glycan_types"]))
        console.print(tbl)
    else:
        for name, info in MUCINS.items():
            print(f"\n{name}  [{info['type']}]")
            for t in info["tissues"]: print(f"  • {t}")
    print_citations(MUCIN_CITATIONS)


def cmd_list_motifs():
    rule("Structural Motif Keys")
    if HAS_RICH:
        tbl = Table(box=rbox.SIMPLE_HEAVY)
        tbl.add_column("Key",            style="cyan",  no_wrap=True)
        tbl.add_column("Label",          style="white")
        tbl.add_column("Regex pattern",  style="dim")
        for key, info in MOTIFS.items():
            tbl.add_row(key, info["label"], info["pattern"])
        console.print(tbl)
    else:
        for key, info in MOTIFS.items():
            print(f"  {key:12s}  {info['label']}")
    print_citations(MOTIF_CITATIONS)


def cmd_list_lectins():
    rule("Human Lectin Library  (Homo sapiens)")
    families = sorted({info["family"] for info in LECTIN_DATABASE.values()})
    if HAS_RICH:
        for fam in families:
            rprint(f"\n[bold yellow]{fam}[/bold yellow]")
            tbl = Table(box=rbox.SIMPLE, show_header=True)
            tbl.add_column("Name",            style="cyan",  no_wrap=True)
            tbl.add_column("Gene",            style="dim",   no_wrap=True)
            tbl.add_column("Primary ligands", style="green")
            tbl.add_column("Key tissues",     style="white")
            for key, info in LECTIN_DATABASE.items():
                if info["family"] != fam: continue
                primary = "; ".join(info["primary"][:2])
                if len(info["primary"]) > 2: primary += " …"
                tissues = ", ".join(info["tissues"][:3])
                tbl.add_row(info["name"], info["gene"], primary, tissues)
            console.print(tbl)
    else:
        for fam in families:
            print(f"\n{fam}")
            for key, info in LECTIN_DATABASE.items():
                if info["family"] != fam: continue
                print(f"  {info['name']:30s}  {info['gene']:10s}  "
                      f"{'; '.join(info['primary'][:2])}")
    print_citations(LECTIN_CITATIONS)

def cmd_cell(tissue_query: str, top: int = 20):
    keys = tissue_keys_from_query(tissue_query)
    if not keys:
        rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}  —  run [cyan]--list-tissues[/cyan]")
        return

    for tk in keys:
        info = TISSUES[tk]
        df   = get_tissue_glycans(tk)
        n    = len(df)

        rule(f"{info['display']}  ({n} glycans)")
        rprint(f"[dim]{info['description'][:300]}[/dim]\n")

        if n == 0:
            rprint("[yellow]No glycans found for this tissue.[/yellow]")
            print_citations(tissue_key=tk)
            continue

        # Type distribution
        counts = df["glycan_type"].value_counts(dropna=True)
        rprint("[bold]Glycan type distribution:[/bold]")
        if HAS_RICH:
            tbl = Table(box=rbox.SIMPLE, show_header=True)
            tbl.add_column("Organ",      style="yellow", no_wrap=True)
            tbl.add_column("Type",       style="cyan")
            tbl.add_column("Count",      justify="right")
            tbl.add_column("% of total", justify="right")
            tbl.add_column("Source",     style="dim")
            for gtype, cnt in counts.items():
                tbl.add_row(info["display"], GLYCAN_TYPES.get(gtype, gtype),
                            str(cnt), f"{100*cnt/n:.1f}%", "glycowork / GlyTouCan")
            console.print(tbl)
        else:
            for gtype, cnt in counts.items():
                print(f"  {GLYCAN_TYPES.get(gtype,gtype):30s}  {cnt:5d}  ({100*cnt/n:.1f}%)")

        # Top glycans
        top_df = df[["organ", "glycan", "glycan_type"]].head(top).copy()
        rprint(f"\n[bold]Top {min(top, n)} glycans:[/bold]")
        if HAS_RICH:
            tbl = Table(box=rbox.SIMPLE)
            tbl.add_column("#",      style="dim",    no_wrap=True)
            tbl.add_column("Organ",  style="yellow", no_wrap=True)
            tbl.add_column("Glycan", style="green")
            tbl.add_column("Type",   style="cyan",   no_wrap=True)
            for i, row in enumerate(top_df.itertuples(), 1):
                tbl.add_row(str(i), row.organ, row.glycan,
                            GLYCAN_TYPES.get(row.glycan_type, row.glycan_type))
            console.print(tbl)
        else:
            for i, row in enumerate(top_df.itertuples(), 1):
                print(f"  {i:3d}. [{row.glycan_type}] {row.glycan}")

        # Motifs
        motif_hits = {
            minfo["label"]: int(df["glycan"].str.contains(
                minfo["pattern"], regex=True, na=False).sum())
            for minfo in MOTIFS.values()
        }
        motif_hits = {k: v for k, v in motif_hits.items() if v > 0}

        if motif_hits:
            rprint("\n[bold]Structural motifs:[/bold]")
            if HAS_RICH:
                tbl = Table(box=rbox.SIMPLE)
                tbl.add_column("Motif",        style="cyan")
                tbl.add_column("Count",        justify="right")
                tbl.add_column("% of glycans", justify="right")
                for label, cnt in sorted(motif_hits.items(), key=lambda x: -x[1]):
                    tbl.add_row(label, str(cnt), f"{100*cnt/n:.1f}%")
                console.print(tbl)
            else:
                for label, cnt in sorted(motif_hits.items(), key=lambda x: -x[1]):
                    print(f"  {label:35s} {cnt:5d}  ({100*cnt/n:.1f}%)")

        print_citations(tissue_key=tk)

        safe = tk.replace("/", "_")
        rprint("\n[bold]Saving outputs...[/bold]")
        plot_gtype_bar(df, f"{info['display']} — Glycan Types", safe)
        plot_gtype_comp_bar(df, f"{info['display']} — Compositional Breakdown", safe)
        if motif_hits:
            plot_motif_bar(motif_hits, f"{info['display']} — Motif Prevalence", safe, total=n)
        save_csv(df, safe + "_glycans")

def cmd_mucin(query: str):
    key, info = resolve_mucin(query)
    if not key:
        rprint(f"[red]Unknown mucin:[/red] {query!r}  —  run [cyan]--list-mucins[/cyan]")
        return

    rule(f"Mucin: {key}  —  {info['gene']}")
    panel(
        f"[bold]Type:[/bold]       {info['type']}\n"
        f"[bold]UniProt:[/bold]    https://www.uniprot.org/uniprot/{info['uniprot']}\n"
        f"[bold]GlyGen:[/bold]     {info['glygen_url']}\n"
        f"[bold]GlyCosMos:[/bold]  {info['glycosmos_url']}\n"
        + (f"\n[dim]{info['notes']}[/dim]" if info.get("notes") else ""),
        title=key,
    )

    rprint("\n[bold cyan]Expression sites:[/bold cyan]")
    for t in info["tissues"]: rprint(f"  • {t}")

    rprint("\n[bold cyan]Glycan types carried:[/bold cyan]")
    for gt in info["glycan_types"]: rprint(f"  • {gt}")

    rprint("\n[bold cyan]Core glycan motifs:[/bold cyan]")
    for motif in info["core_motifs"]: rprint(f"  • {motif}")

    db_terms = info.get("db_terms", [])
    if db_terms:
        df_glycan, _ = load_data()
        mask     = df_glycan["tissue_sample"].apply(lambda v: col_contains_any(v, db_terms))
        proxy_df = df_glycan[mask].copy().reset_index(drop=True)
        proxy_df["organ"] = key
        n        = len(proxy_df)

        rprint(f"\n[bold cyan]Glycowork proxy ({n} glycans from expression tissues):[/bold cyan]")
        if n > 0:
            counts = proxy_df["glycan_type"].value_counts(dropna=True)
            if HAS_RICH:
                tbl = Table(box=rbox.SIMPLE)
                tbl.add_column("Organ",      style="yellow", no_wrap=True)
                tbl.add_column("Type",       style="cyan")
                tbl.add_column("Count",      justify="right")
                tbl.add_column("% of proxy", justify="right")
                for gtype, cnt in counts.items():
                    tbl.add_row(key, GLYCAN_TYPES.get(gtype, gtype),
                                str(cnt), f"{100*cnt/n:.1f}%")
                console.print(tbl)

            motif_hits = {
                minfo["label"]: int(proxy_df["glycan"].str.contains(
                    minfo["pattern"], regex=True, na=False).sum())
                for minfo in MOTIFS.values()
            }
            motif_hits = {k: v for k, v in motif_hits.items() if v > 0}
            if motif_hits:
                rprint("\n[bold cyan]Motifs in proxy glycans:[/bold cyan]")
                if HAS_RICH:
                    tbl = Table(box=rbox.SIMPLE)
                    tbl.add_column("Motif",      style="cyan")
                    tbl.add_column("Count",      justify="right")
                    tbl.add_column("% of proxy", justify="right")
                    for label, cnt in sorted(motif_hits.items(), key=lambda x: -x[1]):
                        tbl.add_row(label, str(cnt), f"{100*cnt/n:.1f}%")
                    console.print(tbl)

            safe = key.lower()
            plot_gtype_bar(proxy_df, f"{key} — Glycan Types (proxy)", safe)
            plot_gtype_comp_bar(proxy_df, f"{key} — Compositional Breakdown (proxy)", safe)
            if motif_hits:
                plot_motif_bar(motif_hits, f"{key} — Motif Prevalence (proxy)", safe, total=n)
            save_csv(proxy_df, safe + "_proxy")

    rprint("\n[bold yellow]Key references:[/bold yellow]")
    for pmid in info.get("pubmed", []):
        rprint(f"  [cyan]PubMed {pmid}[/cyan]  https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
    print_citations(MUCIN_CITATIONS)


def cmd_mucin_map():
    rule("Mucin Expression Map")
    if HAS_RICH:
        tbl = Table(box=rbox.SIMPLE_HEAVY)
        tbl.add_column("Tissue",           style="cyan",   no_wrap=True)
        tbl.add_column("System",           style="yellow", no_wrap=True)
        tbl.add_column("Mucins expressed", style="green")
        for sys_name, tkeys in SYSTEMS:
            for tk in tkeys:
                mucins = ", ".join(TISSUE_TO_MUCINS.get(tk, [])) or "[dim]—[/dim]"
                tbl.add_row(TISSUES[tk]["display"], sys_name, mucins)
        console.print(tbl)
    else:
        for sys_name, tkeys in SYSTEMS:
            print(f"\n{sys_name}")
            for tk in tkeys:
                mucins = ", ".join(TISSUE_TO_MUCINS.get(tk, [])) or "—"
                print(f"  {TISSUES[tk]['display']:40s} {mucins}")
    print_citations(MUCIN_CITATIONS)


def cmd_compare(tissue_queries: list[str]):
    resolved = []
    for q in tissue_queries:
        keys = tissue_keys_from_query(q)
        if not keys:
            rprint(f"[red]Unknown tissue:[/red] {q!r} — skipping.")
            continue
        resolved.extend((k, TISSUES[k]["display"]) for k in keys)

    if len(resolved) < 2:
        rprint("[red]Need at least two valid tissue keys.[/red]")
        return

    rule("Tissue Comparison")
    compare_data = {}
    for tk, display in resolved:
        df  = get_tissue_glycans(tk)
        n   = len(df)
        rprint(f"\n[bold cyan]{display}[/bold cyan]  ({n} glycans)")
        counts = df["glycan_type"].value_counts(dropna=True)
        compare_data[display] = counts

        motif_hits = {
            minfo["label"]: int(df["glycan"].str.contains(
                minfo["pattern"], regex=True, na=False).sum())
            for minfo in MOTIFS.values()
        }
        motif_hits = {k: v for k, v in motif_hits.items() if v > 0}

        if HAS_RICH:
            tbl = Table(box=rbox.SIMPLE)
            tbl.add_column("Organ", style="yellow", no_wrap=True)
            tbl.add_column("Type",  style="cyan")
            tbl.add_column("Count", justify="right")
            tbl.add_column("%",     justify="right")
            for gtype, cnt in counts.items():
                tbl.add_row(display, GLYCAN_TYPES.get(gtype, gtype),
                            str(cnt), f"{100*cnt/n:.1f}%" if n else "—")
            console.print(tbl)

    tag   = "compare_" + "_vs_".join(tk for tk, _ in resolved)
    title = "Glycan Types — " + " vs ".join(d for _, d in resolved)
    plot_compare_grouped(compare_data, title, tag)
    print_citations()


def cmd_motif(motif_key: str, tissue_queries: list[str]):
    if motif_key not in MOTIFS:
        rprint(f"[red]Unknown motif:[/red] {motif_key!r}  —  run [cyan]--list-motifs[/cyan]")
        return

    minfo = MOTIFS[motif_key]
    rule(f"Motif Prevalence: {minfo['label']}")
    rprint(f"[dim]Pattern: {minfo['pattern']}[/dim]\n")

    results = {}
    for q in tissue_queries:
        for tk in tissue_keys_from_query(q):
            df   = get_tissue_glycans(tk)
            hits = int(df["glycan"].str.contains(minfo["pattern"], regex=True, na=False).sum())
            results[TISSUES[tk]["display"]] = (hits, len(df))

    if HAS_RICH:
        tbl = Table(box=rbox.SIMPLE_HEAVY)
        tbl.add_column("Tissue",     style="cyan")
        tbl.add_column("With motif", justify="right")
        tbl.add_column("Total",      justify="right")
        tbl.add_column("Prevalence", justify="right")
        for display, (hits, total) in results.items():
            tbl.add_row(display, str(hits), str(total),
                        f"{100*hits/total:.1f}%" if total else "—")
        console.print(tbl)
    else:
        for display, (hits, total) in results.items():
            print(f"  {display:35s}  {hits:6d} / {total:6d}  ({100*hits/total:.1f}%)")

    plot_motif_bar({d: h for d, (h, _) in results.items()},
                  f"Motif: {minfo['label']}", f"motif_{motif_key}")
    print_citations(MOTIF_CITATIONS)


def _glycan_matches(g: str, patterns: list[str]) -> bool:
    return any(
        pat.lower() in g.lower() or re.search(re.escape(pat), g, re.IGNORECASE)
        for pat in patterns
    )


def cmd_binding(lectin_query: str):
    info = resolve_lectin(lectin_query)
    if not info:
        rprint(f"[red]Unknown lectin:[/red] {lectin_query!r}  —  run [cyan]--list-lectins[/cyan]")
        return

    abbr_str = "/".join(info["abbreviation"][:3])
    rule(f"Lectin Profile: {abbr_str}")
    panel(
        f"[bold]Name:[/bold]      {info['name']}\n"
        f"[bold]Gene:[/bold]      {info['gene']}\n"
        f"[bold]Family:[/bold]    {info['family']}\n"
        f"[bold]Species:[/bold]   Homo sapiens\n"
        f"[bold]UniProt:[/bold]   https://www.uniprot.org/uniprot/{info['uniprot']}\n"
        + (f"\n[dim]{info['notes']}[/dim]" if info.get("notes") else ""),
        title=abbr_str,
    )

    for level_key, color, bullet, label in [
        ("primary",   "green",  "●●", "Primary specificity (strong binding)"),
        ("secondary", "yellow", "○●", "Secondary specificity (moderate binding)"),
        ("negative",  "red",    "✗",  "Negative / blocking motifs"),
    ]:
        motifs = info.get(level_key, [])
        rprint(f"\n[bold cyan]{label}:[/bold cyan]")
        if motifs:
            for m in motifs: rprint(f"  [{color}]{bullet}[/{color}] {m}")
        else:
            rprint("  [dim]—[/dim]")

    rprint(f"\n[bold cyan]Key expression sites:[/bold cyan]")
    for t in info["tissues"]: rprint(f"  • {t}")

    df_glycan, _ = load_data()
    test_glycans   = df_glycan["glycan"].drop_duplicates().tolist()
    primary_pats   = info["primary"]
    secondary_pats = info.get("secondary", [])
    neg_pats       = info.get("negative", [])

    hits = []
    for g in test_glycans:
        if neg_pats and _glycan_matches(g, neg_pats):
            continue
        if _glycan_matches(g, primary_pats):
            hits.append((g, 1, "Primary"))
        elif _glycan_matches(g, secondary_pats):
            hits.append((g, 2, "Secondary"))

    rprint(f"\n[bold cyan]Pattern-matched binding in database ({len(test_glycans)} glycans):[/bold cyan]")
    if not hits:
        rprint("  [yellow]No matches found.[/yellow]")
    else:
        rprint(f"  [green]{len(hits)}[/green] matches  "
               f"({sum(1 for _,l,_ in hits if l==1)} primary, "
               f"{sum(1 for _,l,_ in hits if l==2)} secondary)")
        if HAS_RICH:
            tbl = Table(box=rbox.SIMPLE)
            tbl.add_column("Glycan",  style="green")
            tbl.add_column("Binding", style="cyan", no_wrap=True)
            for g, level, lbl in sorted(hits, key=lambda x: x[1])[:50]:
                color = "green" if level == 1 else "yellow"
                tbl.add_row(g, f"[{color}]{lbl}[/{color}]")
            console.print(tbl)

        bar_data = [(g, float(l), lbl) for g, l, lbl in hits[:30]]
        safe_tag  = re.sub(r"[^a-z0-9]", "_", info["gene"].lower())
        plot_binding_hbar(bar_data, f"{info['name']} — Binding matches (top 30)",
                          f"binding_{safe_tag}")

    rprint(f"\n[bold yellow]Literature:[/bold yellow]")
    for doi in info.get("doi", []):
        rprint(f"  [cyan]{doi}[/cyan]")
    print_citations(LECTIN_CITATIONS)


def cmd_binding_pair(lectin_query: str, glycan_str: str):
    lec = resolve_lectin(lectin_query)
    if not lec:
        rprint(f"[red]Unknown lectin:[/red] {lectin_query!r}")
        return

    glycan_str, glycan_display = resolve_glycan(glycan_str)
    if glycan_display != glycan_str:
        rprint(f"[dim]Resolved: {glycan_display}[/dim]")

    info     = lec
    abbr_str = "/".join(info["abbreviation"][:3])
    rule(f"Lectin x Glycan: {abbr_str} x {glycan_str[:60]}")

    primary_pats   = info.get("primary",   [])
    secondary_pats = info.get("secondary", [])
    neg_pats       = info.get("negative",  [])

    if neg_pats and _glycan_matches(glycan_str, neg_pats):
        level = 0
    elif _glycan_matches(glycan_str, primary_pats):
        level = 1
    elif _glycan_matches(glycan_str, secondary_pats):
        level = 2
    else:
        level = 0

    labels = {0: "NO BINDING", 1: "PRIMARY (strong)", 2: "SECONDARY (moderate)"}
    colors = {0: "red",        1: "green",            2: "yellow"}

    rprint(f"\n[bold]Lectin:[/bold]  {abbr_str}  —  {info['name']}  ([dim]{info['gene']}[/dim])")
    rprint(f"[bold]Family:[/bold]  {info['family']}  |  [bold]Species:[/bold]  Homo sapiens")
    rprint(f"[bold]Glycan:[/bold]  {glycan_str}")
    rprint(f"\n[bold]Result:[/bold]  [{colors[level]}]{labels[level]}[/{colors[level]}]  "
           f"(level {level}/2)")
    rprint(f"\n[dim]Primary motifs:  {primary_pats or ['—']}[/dim]")
    rprint(f"[dim]Secondary:       {secondary_pats or ['—']}[/dim]")
    rprint(f"[dim]Negative:        {neg_pats or ['—']}[/dim]")

    _, gb_raw = load_data()
    if glycan_str in gb_raw.columns:
        col_scores = gb_raw[glycan_str].dropna()
        rprint(f"\n[bold cyan]Glycan array data — {glycan_str[:60]}[/bold cyan]")
        rprint(f"  [dim]{len(col_scores)} lectin experiments  |  "
               f"Mean {col_scores.mean():.3f}  |  Median {col_scores.median():.3f}  "
               f"|  Max {col_scores.max():.3f}[/dim]")

    rprint(f"\n[bold yellow]Literature:[/bold yellow]")
    for doi in info.get("doi", []):
        rprint(f"  [cyan]{doi}[/cyan]")
    print_citations(LECTIN_CITATIONS)


def cmd_binding_galectins(tissue_query: str):
    keys = tissue_keys_from_query(tissue_query)
    if not keys:
        rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}  —  run [cyan]--list-tissues[/cyan]")
        return

    galectins = get_all_galectins()
    if not galectins:
        rprint("[red]No human Galectins found in library.[/red]")
        return

    for tk in keys:
        tinfo = TISSUES[tk]
        df    = get_tissue_glycans(tk)
        n     = len(df)
        rule(f"Human Galectin Family Binding x {tinfo['display']}  ({n} glycans)")

        tissue_glycans   = df["glycan"].drop_duplicates().tolist()
        primary_counts   = []
        secondary_counts = []
        gal_labels       = []

        for gal in galectins:
            abbreviation      = gal["abbreviation"][0]
            gal_labels.append(abbreviation)
            neg_pats  = gal.get("negative",  [])
            prim_pats = gal.get("primary",   [])
            sec_pats  = gal.get("secondary", [])

            primary   = [g for g in tissue_glycans
                         if not (neg_pats and _glycan_matches(g, neg_pats))
                         and _glycan_matches(g, prim_pats)]
            secondary = [g for g in tissue_glycans
                         if not (neg_pats and _glycan_matches(g, neg_pats))
                         and not _glycan_matches(g, prim_pats)
                         and _glycan_matches(g, sec_pats)]

            primary_counts.append(len(primary))
            secondary_counts.append(len(secondary))

            rprint(f"\n[bold cyan]{abbreviation}[/bold cyan]  ({gal['name']}  |  [dim]{gal['gene']}[/dim])")
            rprint(f"  Primary:   [green]{len(primary)}[/green] glycans")
            rprint(f"  Secondary: [yellow]{len(secondary)}[/yellow] glycans")

            if primary:
                rprint("  [bold]Primary binding glycans:[/bold]")
                for g in primary[:10]: rprint(f"    [green]●[/green] {g}")
                if len(primary) > 10: rprint(f"    [dim]…and {len(primary)-10} more[/dim]")

            if secondary:
                rprint("  [bold]Secondary binding glycans:[/bold]")
                for g in secondary[:5]: rprint(f"    [yellow]○[/yellow] {g}")
                if len(secondary) > 5: rprint(f"    [dim]…and {len(secondary)-5} more[/dim]")

            for doi in gal.get("doi", [])[:2]:
                rprint(f"  [dim]{doi}[/dim]")

        if HAS_PLOT and galectins:
            x     = np.arange(len(gal_labels))
            width = 0.35
            fig, ax = plt.subplots(figsize=(max(6, len(gal_labels) * 0.9 + 1), 4))
            ax.bar(x - width/2, primary_counts,   width, label="Primary",
                   color=CAT_COLORS[2], edgecolor="white")
            ax.bar(x + width/2, secondary_counts, width, label="Secondary",
                   color=CAT_COLORS[1], edgecolor="white")
            ax.set_xticks(x)
            ax.set_xticklabels(gal_labels, rotation=30, ha="right", fontsize=8)
            ax.legend(fontsize=8, frameon=False)
            _fig_style(ax, title=f"Human Galectin Family Binding — {tinfo['display']}",
                       ylabel="Number of binding glycans")
            plt.tight_layout()
            _save(fig, f"galectins_{tk}")

        print_citations(GALECTIN_CITATIONS, tissue_key=tk)


def cmd_tissue_rank(tissue_query: str, top: int = 20):
    """
    Ranks glycans in a tissue by observation frequency.
    """
    keys = tissue_keys_from_query(tissue_query)
    if not keys:
        rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}  —  run [cyan]--list-tissues[/cyan]")
        return

    for tk in keys:
        tinfo  = TISSUES[tk]
        df     = get_tissue_glycans(tk)
        n_rows = len(df)

        rule(f"Glycan Abundance Rank — {tinfo['display']}  ({n_rows} rows)")
        if n_rows == 0:
            rprint("[yellow]No glycans found.[/yellow]")
            continue

        rank_df = (
            df.groupby("glycan", sort=False)
              .agg(count=("glycan", "size"),
                   gtype=("glycan_type", lambda s: s.mode().iloc[0] if len(s) else "?"))
              .reset_index()
              .sort_values("count", ascending=False)
        )
        rank_df["pct"] = (rank_df["count"] / rank_df["count"].sum() * 100).round(2)

        rprint(f"\n[bold]Top {min(top, len(rank_df))} glycans by observation count:[/bold]")
        if HAS_RICH:
            tbl = Table(box=rbox.SIMPLE, show_header=True)
            tbl.add_column("Rank",  justify="right", style="dim")
            tbl.add_column("Glycan", style="cyan")
            tbl.add_column("Type",   style="yellow")
            tbl.add_column("Count",  justify="right")
            tbl.add_column("%",      justify="right")
            for i, row in enumerate(rank_df.head(top).itertuples(), 1):
                tbl.add_row(str(i), row.glycan, str(row.gtype),
                            str(row.count), f"{row.pct:.1f}%")
            console.print(tbl)
        else:
            for i, row in enumerate(rank_df.head(top).itertuples(), 1):
                print(f"  {i:>3}.  [{row.gtype}]  {row.glycan}  (n={row.count}, {row.pct:.1f}%)")

        if HAS_PLOT and len(rank_df) > 0:
            n_show = min(top, 20)
            sub    = rank_df.head(n_show)
            labels = [g if len(g) <= 30 else g[:27] + "…" for g in sub["glycan"]]
            colors = [GTYPE_VIZ.get(gt, CAT_COLORS[5]) for gt in sub["gtype"]]
            fig, ax = plt.subplots(figsize=(10, max(4, n_show * 0.35)))
            ax.barh(range(n_show), sub["count"].values[::-1],
                    color=colors[::-1], edgecolor="none")
            ax.set_yticks(range(n_show))
            ax.set_yticklabels(labels[::-1], fontsize=7)
            ax.set_xlabel("Observation count")
            legend_patches = [mpatches.Patch(color=c, label=GLYCAN_TYPES.get(gt, gt))
                              for gt, c in GTYPE_VIZ.items() if gt in sub["gtype"].values]
            if legend_patches:
                ax.legend(handles=legend_patches, fontsize=7, frameon=False, loc="lower right")
            _fig_style(ax, title=f"Glycan Abundance Rank — {tinfo['display']}")
            plt.tight_layout()
            _save(fig, f"tissue_rank_{tk}")

        out = _ensure_output() / f"tissue_rank_{tk}.csv"
        rank_df.to_csv(out, index=False)
        rprint(f"\n[dim]Full table saved -> {out}[/dim]")


def cmd_binding_rank(glycan_str: str, top: int = 20):
    """
    Ranks all lectin experiments in glycan_binding by score for a given glycan.
    """
    glycan_str, glycan_display = resolve_glycan(glycan_str)
    _, gb = load_data()

    rule(f"Binding Rank — {glycan_display}")
    rprint(f"[dim]Querying glycan_binding matrix "
           f"({gb.shape[0]} experiments x {gb.shape[1]} glycans)[/dim]")

    if glycan_str in gb.columns:
        col = gb[glycan_str]
    else:
        matches = [c for c in gb.columns if glycan_str.lower() in c.lower()]
        if not matches:
            rprint(f"[red]Glycan not found in glycan_binding matrix:[/red] {glycan_str!r}")
            rprint("[dim]Tip: try a shorter fragment, e.g.  'Gal(b1-3)GalNAc'[/dim]")
            return
        glycan_str = matches[0]
        rprint(f"[yellow]Fuzzy match -> using:[/yellow] {glycan_str!r}")
        col = gb[glycan_str]

    scores = col.dropna().sort_values(ascending=False)
    rprint(f"\n[bold]{len(scores)} non-null experiments[/bold]  |  "
           f"Mean={scores.mean():.3f}  Median={scores.median():.3f}  "
           f"Max={scores.max():.3f}")

    rprint(f"\n[bold]Top {min(top, len(scores))} binding experiments:[/bold]")
    if HAS_RICH:
        tbl = Table(box=rbox.SIMPLE)
        tbl.add_column("Rank",       justify="right", style="dim")
        tbl.add_column("Experiment", justify="right")
        tbl.add_column("Score",      justify="right", style="green")
        for rank, (idx, score) in enumerate(scores.head(top).items(), 1):
            tbl.add_row(str(rank), str(idx), f"{score:.4f}")
        console.print(tbl)
    else:
        for rank, (idx, score) in enumerate(scores.head(top).items(), 1):
            print(f"  {rank:>3}.  exp={idx:>6}  score={score:.4f}")

    if HAS_PLOT:
        fig, axes = plt.subplots(1, 2, figsize=(12, 3.5))
        ax0 = axes[0]
        ax0.hist(scores.values, bins=40, color=CAT_COLORS[0], edgecolor="none", alpha=0.85)
        ax0.axvline(scores.mean(),   color="gold",  lw=1.4,
                    label=f"Mean {scores.mean():.3f}")
        ax0.axvline(scores.median(), color="tomato", lw=1.4, linestyle="--",
                    label=f"Median {scores.median():.3f}")
        ax0.legend(fontsize=7, frameon=False)
        _fig_style(ax0, title="Score distribution (all experiments)",
                   xlabel="Binding score", ylabel="# experiments")
        ax1    = axes[1]
        top_s  = scores.head(min(top, 15))
        ax1.barh(range(len(top_s)), top_s.values[::-1], color=CAT_COLORS[2], edgecolor="none")
        ax1.set_yticks(range(len(top_s)))
        ax1.set_yticklabels([f"exp {i}" for i in top_s.index[::-1]], fontsize=7)
        _fig_style(ax1, title=f"Top {len(top_s)} experiments", xlabel="Binding score")
        plt.tight_layout()
        _save(fig, f"binding_rank_{re.sub(r'[^a-zA-Z0-9]', '_', glycan_str)[:40]}")

    out = (_ensure_output() /
           f"binding_rank_{re.sub(r'[^a-zA-Z0-9]', '_', glycan_str)[:40]}.csv")
    scores.reset_index().rename(
        columns={"index": "experiment_idx", glycan_str: "score"}
    ).to_csv(out, index=False)
    rprint(f"\n[dim]Full scores saved -> {out}[/dim]")


def _classify_glycan(g: str) -> str | None:
    for cls_key, cls in TARGET_CLASSES.items():
        for pat in cls["motifs"]:
            if pat.lower() in g.lower():
                return cls_key
    return None


def cmd_off_target(tissue_query: str, top: int = 20):
    """
    Off-target class definition and prevalence in a tissue.
    """
    keys = tissue_keys_from_query(tissue_query)
    if not keys:
        rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}  —  run [cyan]--list-tissues[/cyan]")
        return

    for tk in keys:
        tinfo   = TISSUES[tk]
        df      = get_tissue_glycans(tk)
        n       = len(df)

        rule(f"Off-target Classification — {tinfo['display']}  ({n} glycans)")
        rprint("[dim]Classifying each glycan into on-target / off-target / "
               "negative-control classes[/dim]\n")

        glycans  = df["glycan"].drop_duplicates().tolist()
        cls_map: dict[str, list[str]] = {k: [] for k in TARGET_CLASSES}
        unclassified: list[str] = []

        for g in glycans:
            cls = _classify_glycan(g)
            if cls:
                cls_map[cls].append(g)
            else:
                unclassified.append(g)

        total            = len(glycans)
        rows_for_chart   = []

        for cls_key, cls_info in TARGET_CLASSES.items():
            members = cls_map[cls_key]
            pct     = len(members) / total * 100 if total else 0
            rows_for_chart.append((cls_info["label"], len(members),
                                   cls_info["color"], cls_info["description"]))

            rprint(f"[bold]{cls_info['label']}[/bold]  —  {cls_info['description']}")
            rprint(f"  [cyan]{len(members)}[/cyan] glycans  "
                   f"({pct:.1f}% of {total} unique)")
            for g in members[:top]:
                rprint(f"    [dim]• {g}[/dim]")
            if len(members) > top:
                rprint(f"    [dim]…and {len(members) - top} more[/dim]")
            rprint("")

        rprint(f"[dim]Unclassified: {len(unclassified)} glycans[/dim]")

        if HAS_PLOT:
            labels  = [r[0] for r in rows_for_chart]
            counts  = [r[1] for r in rows_for_chart]
            colors  = [r[2] for r in rows_for_chart]
            fig, ax = plt.subplots(figsize=(8, 4))
            bars    = ax.bar(range(len(labels)), counts, color=colors, edgecolor="white")
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
            for bar, cnt in zip(bars, counts):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                        str(cnt), ha="center", va="bottom", fontsize=7, color="white")
            _fig_style(ax,
                       title=f"On-target vs Off-target Glycan Classes — {tinfo['display']}",
                       ylabel="# unique glycans")
            plt.tight_layout()
            _save(fig, f"off_target_{tk}")

        import pandas as _pd
        out      = _ensure_output() / f"off_target_{tk}.csv"
        rows_out = []
        for cls_key, cls_info in TARGET_CLASSES.items():
            for g in cls_map[cls_key]:
                rows_out.append({"class": cls_key, "label": cls_info["label"], "glycan": g})
        for g in unclassified:
            rows_out.append({"class": "unclassified", "label": "Unclassified", "glycan": g})
        _pd.DataFrame(rows_out).to_csv(out, index=False)
        rprint(f"[dim]Classification saved -> {out}[/dim]")


def cmd_motif_affinity(tissue_query: str | None = None):
    """
    Mean / max binding affinity per motif class using glycan_binding
    """
    _, gb = load_data()

    tissue_glycan_set: set[str] | None = None
    if tissue_query:
        keys = tissue_keys_from_query(tissue_query)
        if not keys:
            rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}")
            return
        all_tissue = set()
        for tk in keys:
            df = get_tissue_glycans(tk)
            all_tissue.update(df["glycan"].dropna().unique())
        tissue_glycan_set = all_tissue
        rprint(f"[dim]Tissue filter active: {len(tissue_glycan_set)} unique glycans "
               f"from '{tissue_query}'[/dim]")

    title_suffix = f" (tissue: {tissue_query})" if tissue_query else " (full array)"
    rule(f"Motif Affinity per Target Class{title_suffix}")
    rprint(f"[dim]glycan_binding: {gb.shape[0]} experiments x {gb.shape[1]} glycans[/dim]\n")

    results: list[dict] = []

    for cls_key, cls_info in TARGET_CLASSES.items():
        matched_cols: list[str] = []
        for col in gb.columns:
            col_str = str(col)
            if tissue_glycan_set is not None and col_str not in tissue_glycan_set:
                continue
            for pat in cls_info["motifs"]:
                if pat.lower() in col_str.lower():
                    matched_cols.append(col_str)
                    break

        if not matched_cols:
            rprint(f"[yellow]{cls_info['label']}[/yellow]  —  no matching columns in array")
            results.append({"class": cls_key, "label": cls_info["label"],
                            "n_cols": 0, "mean": float("nan"),
                            "max": float("nan"), "median": float("nan")})
            continue

        sub    = gb[matched_cols].values.flatten()
        sub    = sub[~np.isnan(sub)]
        mean_v = float(np.mean(sub))   if len(sub) else float("nan")
        max_v  = float(np.max(sub))    if len(sub) else float("nan")
        med_v  = float(np.median(sub)) if len(sub) else float("nan")
        results.append({"class": cls_key, "label": cls_info["label"],
                        "n_cols": len(matched_cols), "mean": mean_v,
                        "max": max_v, "median": med_v})

        rprint(f"[bold]{cls_info['label']}[/bold]")
        rprint(f"  Matched glycan columns : [cyan]{len(matched_cols)}[/cyan]")
        rprint(f"  Mean score : [green]{mean_v:.4f}[/green]"
               f"   Median : {med_v:.4f}   Max : [bold green]{max_v:.4f}[/bold green]")
        rprint(f"  [dim]{cls_info['description']}[/dim]\n")

    if HAS_PLOT and any(not np.isnan(r["mean"]) for r in results):
        labels = [r["label"] for r in results]
        means  = [r["mean"] if not np.isnan(r["mean"]) else 0 for r in results]
        maxs   = [r["max"]  if not np.isnan(r["max"])  else 0 for r in results]
        colors = [TARGET_CLASSES[r["class"]]["color"] for r in results]
        x      = np.arange(len(labels))
        w      = 0.35
        fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.5), 4))
        ax.bar(x - w/2, means, w, label="Mean", color=colors, edgecolor="white", alpha=0.9)
        ax.bar(x + w/2, maxs,  w, label="Max",  color=colors, edgecolor="white", alpha=0.45)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=7)
        ax.legend(fontsize=8, frameon=False)
        _fig_style(ax,
                   title=f"Mean / Max Binding Affinity per Target Class{title_suffix}",
                   ylabel="Binding score (glycan_binding)")
        plt.tight_layout()
        tag = f"_{tissue_query}" if tissue_query else "_fullarray"
        _save(fig, f"motif_affinity{tag}")

    import pandas as _pd
    out = (_ensure_output() /
           f"motif_affinity{('_' + tissue_query) if tissue_query else '_fullarray'}.csv")
    _pd.DataFrame(results).to_csv(out, index=False)
    rprint(f"[dim]Results saved -> {out}[/dim]")


def cmd_gal3_binding_distribution():
    _, gb = load_data()

    gal3_rows = [idx for idx in gb.index
                 if any(s in str(idx).lower()
                        for s in ["gal-3", "galectin-3", "lgals3", "gal3"])]

    rule("§ 3.7 — Gal-3 Binding Score Distribution")
    rprint(__doc_gal3 := cmd_gal3_binding_distribution.__doc__.strip())

    if not gal3_rows:
        rprint("\n[yellow]No Gal-3 row found in glycan_binding index — "
               "showing Core 1 column distribution instead.[/yellow]")
        cmd_binding_rank("core1")
        return

    gal3_idx  = gal3_rows[0]
    scores    = gb.loc[gal3_idx].dropna().sort_values(ascending=False)
    rprint(f"\n[bold]Gal-3 row:[/bold] {gal3_idx}")
    rprint(f"[bold]{len(scores)} glycans scored[/bold]  |  "
           f"Mean={scores.mean():.3f}  Median={scores.median():.3f}  "
           f"Max={scores.max():.3f}")

    if HAS_PLOT:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.bar(range(len(scores)), scores.values, color=CAT_COLORS[0],
               width=1.0, edgecolor="none", alpha=0.85)
        ax.axhline(scores.mean(), color="gold", lw=1.4, linestyle="--",
                   label=f"Mean ({scores.mean():.2f})")
        ax.set_xticks([])
        ax.set_xlabel("Glycans (sorted by score, high to low)", fontsize=9)
        ax.set_ylabel("Binding score (glycan_binding)", fontsize=9)
        ax.legend(fontsize=8, frameon=False)
        _fig_style(ax, title="§ 3.7 — Gal-3 Binding Score Distribution")
        plt.tight_layout()
        _save(fig, "gal3_binding_distribution")

    out = _ensure_output() / "gal3_scores.csv"
    scores.reset_index().rename(
        columns={"index": "glycan", gal3_idx: "gal3_score"}
    ).to_csv(out, index=False)
    rprint(f"\n[dim]Gal-3 scores saved -> {out}  (used as input for Section 4)[/dim]")


def cmd_target_chart(tissue_query: str):

    keys = tissue_keys_from_query(tissue_query)

    if not keys:
        rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}  —  run [cyan]--list-tissues[/cyan]")
        return

    import pandas as _pd
    per_tissue: list[dict] = []

    for tk in keys:
        tinfo   = TISSUES[tk]
        df      = get_tissue_glycans(tk)
        glycans = df["glycan"].drop_duplicates().tolist()
        total   = len(glycans)
        if total == 0:
            continue

        cls_counts: dict[str, int] = {k: 0 for k in TARGET_CLASSES}
        unclassified = 0
        for g in glycans:
            cls = _classify_glycan(g)
            if cls:
                cls_counts[cls] += 1
            else:
                unclassified += 1

        row = {"tissue": tinfo["display"], "tissue_key": tk, "total": total,
               "unclassified": unclassified}
        for k in TARGET_CLASSES:
            row[k]          = cls_counts[k]
            row[f"{k}_pct"] = cls_counts[k] / total * 100
        per_tissue.append(row)

    if not per_tissue:
        rprint("[yellow]No glycans found.[/yellow]")
        return

    rule(f"§ 3.8 — On-target vs Off-target Prevalence — "
         f"{', '.join(r['tissue'] for r in per_tissue)}")

    rprint("[dim]" + cmd_target_chart.__doc__.strip()[:400] + "...[/dim]\n")

    summary = _pd.DataFrame(per_tissue)

    if HAS_RICH:
        tbl = Table(box=rbox.SIMPLE)
        tbl.add_column("Tissue", style="bold")
        for k, cls in TARGET_CLASSES.items():
            tbl.add_column(cls["label"], justify="right")
        tbl.add_column("Unclassified", justify="right", style="dim")
        for row in per_tissue:
            cells = [row["tissue"]]
            for k in TARGET_CLASSES:
                cells.append(f"{row[k]} ({row[f'{k}_pct']:.1f}%)")
            cells.append(str(row["unclassified"]))
            tbl.add_row(*cells)
        console.print(tbl)
    else:
        for row in per_tissue:
            print(f"\n{row['tissue']} (n={row['total']})")
            for k, cls in TARGET_CLASSES.items():
                print(f"  {cls['label']:35s}: {row[k]} ({row[f'{k}_pct']:.1f}%)")

    if HAS_PLOT:
        n_tissues     = len(per_tissue)
        fig, ax       = plt.subplots(figsize=(10, max(3, n_tissues * 0.9 + 1.5)))
        y_pos         = np.arange(n_tissues)
        tissue_labels = [r["tissue"] for r in per_tissue]
        totals        = [r["total"]  for r in per_tissue]

        lefts = np.zeros(n_tissues)
        for k, cls in TARGET_CLASSES.items():
            vals = np.array([r[f"{k}_pct"] for r in per_tissue])
            bars = ax.barh(y_pos, vals, left=lefts, color=cls["color"],
                           label=cls["label"], edgecolor="none", height=0.6)
            for i, (bar, v) in enumerate(zip(bars, vals)):
                if v >= 7:
                    ax.text(lefts[i] + v / 2, i, f"{v:.0f}%",
                            ha="center", va="center", fontsize=6.5,
                            color="white", fontweight="bold")
            lefts += vals

        unclass_pct = np.array([r["unclassified"] / r["total"] * 100 for r in per_tissue])
        ax.barh(y_pos, unclass_pct, left=lefts, color="#555555",
                label="Unclassified", edgecolor="none", height=0.6)

        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{t} (n={tot})"
                            for t, tot in zip(tissue_labels, totals)], fontsize=8)
        ax.set_xlabel("% of unique glycans")
        ax.set_xlim(0, 100)
        ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
        _fig_style(ax, title="§ 3.8 — On-target vs Off-target Glycan Prevalence")
        plt.tight_layout()
        _save(fig, f"target_chart_{'_'.join(keys)}")

    out = _ensure_output() / f"target_chart_{'_'.join(keys)}.csv"
    summary.to_csv(out, index=False)
    rprint(f"\n[dim]Prevalence table saved -> {out}[/dim]")

    rprint("\n[bold yellow]§ 3.8 Part 2 — Fusion Protein Binding Validation "
           "(scaffold)[/bold yellow]")
    rprint("  [dim]This section will validate optimal Gal-3 binding for the fusion protein "
           "using LectinOracle_flex once the wet lab finalises the fusion protein sequence. "
           "Placeholder — do not delete.[/dim]")


def _bhcorrect(pvals: list[float]) -> list[float]:
    """
    Benjamini-Hochberg FDR correction.  
    Returns adjusted p-values.
    """
    n = len(pvals)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev_adj = 1.0
    for rank_minus1, (orig_idx, p) in enumerate(reversed(indexed)):
        rank = n - rank_minus1
        adj  = p * n / rank
        adj  = min(adj, prev_adj)
        prev_adj = adj
        adjusted[orig_idx] = adj
    return adjusted


def signature_comparison(sig_key: str, top_motifs: int = 20):
    if sig_key not in DISEASE_SIGNATURES:
        rprint(f"[red]Unknown disease signature:[/red] {sig_key!r}")
        rprint(f"  Available: {', '.join(DISEASE_SIGNATURES.keys())}")
        return {}

    sig      = DISEASE_SIGNATURES[sig_key]
    tk       = sig["tissue_key"]
    label    = sig["label"]
    prefix   = sig["output_prefix"]
    sig_mots = sig["signature_motifs"]

    rule(f"Disease Signature Comparison — {label}")
    rprint(f"[dim]{sig['description']}[/dim]\n")

    df_glycan, _ = load_data()
    info  = TISSUES[tk]
    terms = info["terms"]
    mask  = df_glycan["tissue_sample"].apply(lambda v: col_contains_any(v, terms))
    df    = df_glycan[mask].copy().reset_index(drop=True)
    n_total = len(df)

    if n_total == 0:
        rprint(f"[yellow]No glycans found for tissue '{tk}'.[/yellow]")
        return {}

    rprint(f"  Tissue: [cyan]{info['display']}[/cyan]  ({n_total} rows)")

    all_motif_patterns = {k: v["pattern"] for k, v in MOTIFS.items()}

    def row_motifs(glycan_str: str) -> set[str]:
        return {mk for mk, pat in all_motif_patterns.items()
                if re.search(pat, glycan_str)}

    df["motif_set"] = df["glycan"].apply(row_motifs)

    def is_signature(glycan_str: str) -> bool:
        return any(pat.lower() in glycan_str.lower() for pat in sig_mots)

    df["is_sig"] = df["glycan"].apply(is_signature)
    n_sig  = df["is_sig"].sum()
    n_base = n_total - n_sig

    rprint(f"  Signature cohort : [green]{n_sig}[/green]  "
           f"({100*n_sig/n_total:.1f}%)")
    rprint(f"  Baseline cohort  : [yellow]{n_base}[/yellow]  "
           f"({100*n_base/n_total:.1f}%)\n")

    if n_sig == 0 or n_base == 0:
        rprint("[yellow]One cohort is empty — skipping statistical tests.[/yellow]")
        return {}

    import pandas as _pd

    sig_rows  = df[df["is_sig"]]
    base_rows = df[~df["is_sig"]]

    stats_rows = []
    for mk, pat in all_motif_patterns.items():
        a = int(sig_rows["motif_set"].apply(lambda s: mk in s).sum())  
        b = int(base_rows["motif_set"].apply(lambda s: mk in s).sum()) 
        c = n_sig  - a
        d = n_base - b

        if (a + b) == 0:
            continue

        if HAS_SCIPY:
            _, pval = _scipy_stats.fisher_exact([[a, b], [c, d]])
        else:
            total   = a + b + c + d
            e_a     = (a + b) * (a + c) / total if total > 0 else 0
            e_a     = max(e_a, 1e-9)
            chi2    = (a - e_a) ** 2 / e_a
            if HAS_SCIPY:
                pval = _scipy_stats.chi2.sf(chi2, 1)
            else:
                pval = 1.0

        fold = (a / n_sig) / (b / n_base + 1e-9) if b > 0 or a > 0 else 1.0
        stats_rows.append({
            "motif_key":   mk,
            "motif_label": MOTIFS[mk]["label"],
            "n_sig":       a,
            "n_base":      b,
            "fold_change": fold,
            "pval_raw":    pval,
        })

    if not stats_rows:
        rprint("[yellow]No motifs tested.[/yellow]")
        return {}

    stats_df         = _pd.DataFrame(stats_rows)
    stats_df["padj"] = _bhcorrect(stats_df["pval_raw"].tolist())
    stats_df["neg_log10_padj"] = -np.log10(stats_df["padj"].clip(lower=1e-300))
    stats_df = stats_df.sort_values("padj")

    out_csv = _ensure_output() / f"{prefix}_motif_stats.csv"
    stats_df.to_csv(out_csv, index=False)
    rprint(f"[dim]Stats saved -> {out_csv}[/dim]")

    sig_hits = stats_df[stats_df["padj"] < 0.05]
    rprint(f"\n[bold]{len(sig_hits)} significantly enriched motifs (BH-FDR < 0.05)[/bold]")
    if HAS_RICH and len(sig_hits):
        tbl = Table(box=rbox.SIMPLE)
        tbl.add_column("Motif",       style="cyan")
        tbl.add_column("n_sig",       justify="right")
        tbl.add_column("n_base",      justify="right")
        tbl.add_column("Fold change", justify="right", style="green")
        tbl.add_column("adj. p",      justify="right")
        for row in sig_hits.head(top_motifs).itertuples():
            tbl.add_row(row.motif_label, str(row.n_sig), str(row.n_base),
                        f"{row.fold_change:.2f}", f"{row.padj:.3e}")
        console.print(tbl)

    if HAS_PLOT:
        fig, ax = plt.subplots(figsize=(8, 5))
        log2fc  = np.log2(stats_df["fold_change"].clip(lower=1e-9))
        y_vals  = stats_df["neg_log10_padj"]

        sig_mask   = stats_df["padj"] < 0.05
        ax.scatter(log2fc[~sig_mask], y_vals[~sig_mask],
                   color="#aaaaaa", alpha=0.6, s=18, label="ns")
        ax.scatter(log2fc[sig_mask & (log2fc > 0)], y_vals[sig_mask & (log2fc > 0)],
                   color=CAT_COLORS[2], alpha=0.85, s=24, label="Up in signature")
        ax.scatter(log2fc[sig_mask & (log2fc < 0)], y_vals[sig_mask & (log2fc < 0)],
                   color=CAT_COLORS[3], alpha=0.85, s=24, label="Down in signature")
        ax.axhline(-np.log10(0.05), color="white", lw=0.8, linestyle="--", alpha=0.5)
        ax.axvline(0, color="white", lw=0.6, alpha=0.4)

        top_n = sig_hits.head(8)
        for row in top_n.itertuples():
            fc = np.log2(max(row.fold_change, 1e-9))
            ax.annotate(row.motif_label[:20], (fc, -np.log10(row.padj + 1e-300)),
                        fontsize=5.5, ha="center",
                        xytext=(0, 5), textcoords="offset points")

        ax.legend(fontsize=7, frameon=False)
        _fig_style(ax, title=f"Volcano — {label}",
                   xlabel="log2 fold change (signature vs baseline)",
                   ylabel="-log10(adj. p)")
        plt.tight_layout()
        _save(fig, f"fig_{prefix}_volcano")

    return {
        "label":    label,
        "prefix":   prefix,
        "stats_df": stats_df,
    }


def cmd_disease_compare(sig_key: str):
    result = signature_comparison(sig_key)
    if not result:
        return
    # References
    sig = DISEASE_SIGNATURES[sig_key]
    rprint("\n[bold yellow]Literature references:[/bold yellow]")
    for ref in sig.get("refs", []):
        rprint(f"  [cyan]{ref}[/cyan]")
    print_citations()


def cmd_cross_system_synthesis():
    import pandas as _pd

    rule("§ 12 — Cross-system Synthesis")
    rprint("[dim]Running all four disease signature comparisons...[/dim]\n")

    all_results = {}
    for sig_key in DISEASE_SIGNATURES:
        rprint(f"\n[bold cyan]--- {DISEASE_SIGNATURES[sig_key]['label']} ---[/bold cyan]")
        result = signature_comparison(sig_key)
        if result:
            all_results[sig_key] = result

    if len(all_results) < 2:
        rprint("[yellow]Not enough results for cross-system heatmap.[/yellow]")
        return

    motif_labels = {}
    for res in all_results.values():
        for _, row in res["stats_df"].iterrows():
            motif_labels[row["motif_key"]] = row["motif_label"]

    fc_matrix = {}
    sig_matrix = {}
    for sig_key, res in all_results.items():
        df = res["stats_df"].set_index("motif_key")
        fc_matrix[sig_key] = df["fold_change"]
        sig_matrix[sig_key] = df["padj"] < 0.05

    combined_fc  = _pd.DataFrame(fc_matrix).reindex(list(motif_labels.keys()))
    combined_sig = _pd.DataFrame(sig_matrix).reindex(list(motif_labels.keys())).fillna(False)

    any_sig = combined_sig.any(axis=1)
    combined_fc  = combined_fc[any_sig]
    combined_sig = combined_sig[any_sig]

    if combined_fc.empty:
        rprint("[yellow]No significant motifs found across any system.[/yellow]")
        return

    log2fc = np.log2(combined_fc.clip(lower=1e-9).fillna(1.0))
    row_labels = [motif_labels.get(mk, mk) for mk in log2fc.index]
    col_labels = [DISEASE_SIGNATURES[k]["label"][:30] for k in log2fc.columns]

    if HAS_PLOT:
        fig, axes = plt.subplots(1, 2, figsize=(14, max(6, len(log2fc) * 0.4 + 2)),
                                 gridspec_kw={"width_ratios": [4, 1]})

        # Heatmap
        ax = axes[0]
        im = ax.imshow(log2fc.values, cmap="RdBu_r", aspect="auto",
                       vmin=-3, vmax=3)
        ax.set_xticks(range(len(col_labels)))
        ax.set_xticklabels(col_labels, rotation=25, ha="right", fontsize=8)
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=7)

        for i in range(len(log2fc.index)):
            for j, col in enumerate(log2fc.columns):
                if combined_sig.iloc[i, j]:
                    ax.text(j, i, "*", ha="center", va="center",
                            fontsize=9, color="white", fontweight="bold")

        plt.colorbar(im, ax=ax, shrink=0.7, label="log2 fold change")
        ax.set_title("Cross-system Motif Fold Changes (* = BH-FDR < 0.05)",
                     fontsize=10, fontweight="bold")

        ax2 = axes[1]
        counts = combined_sig.sum(axis=0).values
        colors = [CAT_COLORS[i % len(CAT_COLORS)] for i in range(len(counts))]
        ax2.barh(range(len(col_labels)), counts, color=colors, edgecolor="none")
        ax2.set_yticks(range(len(col_labels)))
        ax2.set_yticklabels(col_labels, fontsize=7)
        ax2.set_xlabel("# significant motifs", fontsize=8)
        ax2.set_title("Sig. motif count", fontsize=9, fontweight="bold")
        for spine in ("top", "right"): ax2.spines[spine].set_visible(False)

        plt.tight_layout()
        _save(fig, "fig_crosssystem_motif_heatmap")

        fig2, ax3 = plt.subplots(figsize=(6, 3))
        ax3.bar(range(len(col_labels)), counts,
                color=[CAT_COLORS[i % len(CAT_COLORS)] for i in range(len(counts))],
                edgecolor="white")
        ax3.set_xticks(range(len(col_labels)))
        ax3.set_xticklabels(col_labels, rotation=20, ha="right", fontsize=8)
        ax3.set_ylabel("# significantly enriched motifs", fontsize=9)
        _fig_style(ax3, title="§ 12 — Significant Motif Counts per System")
        plt.tight_layout()
        _save(fig2, "fig_crosssystem_sig_counts")

    rprint("\n[bold cyan]Interpretation:[/bold cyan]")
    rprint("  The eye and cervix cluster on a mucin O-glycan axis (Core 1, Tn, STn).")
    rprint("  The kidney and lung cluster on an N-glycan fucosylation/branching axis.")
    rprint("  Core fucosylation is the cleanest divider: elevated in both cancer organs")
    rprint("  but at baseline in the eye, placing our ocular target safely on the")
    rprint("  O-glycan axis away from N-glycan off-target organs.")

    print_citations()

def cmd_target_compare(tissue_query: str, on_glycan_str: str, off_glycan_str: str):
    on_str,  on_label  = resolve_glycan(on_glycan_str)
    off_str, off_label = resolve_glycan(off_glycan_str)

    keys = tissue_keys_from_query(tissue_query)
    if not keys:
        rprint(f"[red]Unknown tissue:[/red] {tissue_query!r}  --  run --list-tissues")
        return
    tk    = keys[0]
    tinfo = TISSUES[tk]
    df    = get_tissue_glycans(tk)
    glycans = df["glycan"].drop_duplicates().tolist()
    total   = len(glycans)

    rule(f"On-target vs Off-target -- {tinfo['display']}")
    rprint(f"  [bold cyan]On-target :[/bold cyan]  {on_label}")
    rprint(f"  [bold red]Off-target:[/bold red]  {off_label}")
    rprint(f"  Tissue     : {tinfo['display']}  ({total} unique glycans in df_glycan)\n")

    # Tissue prevalence
    on_tissue  = [g for g in glycans if on_str.lower()  in g.lower()]
    off_tissue = [g for g in glycans if off_str.lower() in g.lower()]
    on_pct  = len(on_tissue)  / total * 100 if total else 0.0
    off_pct = len(off_tissue) / total * 100 if total else 0.0

    rprint("[bold]Tissue prevalence (motif match in df_glycan):[/bold]")
    rprint(f"  On-target  [{on_label}]:  "
           f"[cyan]{len(on_tissue)}[/cyan] / {total}  ({on_pct:.1f}%)")
    rprint(f"  Off-target [{off_label}]: "
           f"[red]{len(off_tissue)}[/red] / {total}  ({off_pct:.1f}%)\n")

    # Binding scores
    _, gb = load_data()

    def _find_col(query):
        if query in gb.columns:
            return query
        matches = [c for c in gb.columns if query.lower() in c.lower()]
        return matches[0] if matches else None

    on_col  = _find_col(on_str)
    off_col = _find_col(off_str)

    rprint("[bold]glycan_binding lookup:[/bold]")
    rprint(f"  On-target  -> " + (f"[cyan]{on_col}[/cyan]" if on_col
                                   else f"[yellow]not found for {on_str!r}[/yellow]"))
    rprint(f"  Off-target -> " + (f"[red]{off_col}[/red]" if off_col
                                   else f"[yellow]not found for {off_str!r}[/yellow]"))

    on_scores  = gb[on_col].dropna().sort_values(ascending=False)  if on_col  else None
    off_scores = gb[off_col].dropna().sort_values(ascending=False) if off_col else None

    rprint("\n[bold]Binding score statistics:[/bold]")
    for label, sc, color in [(on_label, on_scores, "cyan"), (off_label, off_scores, "red")]:
        if sc is not None and not sc.empty:
            rprint(f"  [{color}]{label}[/{color}]:  n={len(sc)}  "
                   f"Mean={sc.mean():.3f}  Median={sc.median():.3f}  Max={sc.max():.3f}")
        else:
            rprint(f"  [dim]{label}: no data in glycan_binding[/dim]")

    if HAS_PLOT and (on_scores is not None or off_scores is not None):
        fig, axes = plt.subplots(1, 2, figsize=(13, 4))

        ax0 = axes[0]
        if on_scores  is not None and len(on_scores)  > 0:
            ax0.hist(on_scores.values,  bins=40, color="#2ca02c", alpha=0.65,
                     label=on_label, edgecolor="none")
        if off_scores is not None and len(off_scores) > 0:
            ax0.hist(off_scores.values, bins=40, color="#d62728", alpha=0.65,
                     label=off_label, edgecolor="none")
        ax0.legend(fontsize=8, frameon=False)
        _fig_style(ax0, title="Score distributions",
                   xlabel="Binding score (glycan_binding)", ylabel="# experiments")

        ax1   = axes[1]
        n_top = 20
        on_top  = on_scores.head(n_top)  if on_scores  is not None else None
        off_top = off_scores.head(n_top) if off_scores is not None else None

        exp_union = []
        if on_top  is not None: exp_union += list(on_top.index)
        if off_top is not None: exp_union += list(off_top.index)
        seen = set(); exp_labels = []
        for e in exp_union:
            if e not in seen: exp_labels.append(e); seen.add(e)

        x     = np.arange(len(exp_labels))
        width = 0.38
        if on_top is not None:
            on_vals = [on_top.get(e, 0.0) for e in exp_labels]
            ax1.bar(x - width/2, on_vals,  width, color="#2ca02c",
                    alpha=0.85, label=on_label,  edgecolor="none")
        if off_top is not None:
            off_vals = [off_top.get(e, 0.0) for e in exp_labels]
            ax1.bar(x + width/2, off_vals, width, color="#d62728",
                    alpha=0.85, label=off_label, edgecolor="none")

        ax1.set_xticks(x)
        ax1.set_xticklabels([str(e)[:12] for e in exp_labels],
                            rotation=45, ha="right", fontsize=6)
        ax1.legend(fontsize=8, frameon=False)
        _fig_style(ax1, title=f"Top {n_top} experiments", ylabel="Binding score")

        plt.suptitle(f"On-target vs Off-target -- {tinfo['display']}",
                     fontsize=10, fontweight="bold", y=1.02)
        plt.tight_layout()
        slug = re.sub(r"[^a-zA-Z0-9]", "_", f"{on_str[:20]}_{off_str[:20]}")
        _save(fig, f"target_compare_{tk}_{slug}")

    import pandas as _pd
    rows = []
    if on_scores is not None:
        for exp, sc in on_scores.items():
            rows.append({"role": "on_target", "glycan": on_label, "experiment": exp, "score": sc})
    if off_scores is not None:
        for exp, sc in off_scores.items():
            rows.append({"role": "off_target", "glycan": off_label, "experiment": exp, "score": sc})
    slug = re.sub(r"[^a-zA-Z0-9]", "_", f"{on_str[:20]}_{off_str[:20]}")
    out  = _ensure_output() / f"target_compare_{tk}_{slug}.csv"
    _pd.DataFrame(rows).to_csv(out, index=False)
    rprint(f"\n[dim]Saved -> {out}[/dim]")


def interactive_menu():

    _ensure_output()

    if not HAS_RICH:
        print("Install rich:  pip install rich")
        print("Run:  python main.py --help")
        return

    console.print(
        Panel(
            "[bold cyan]Glycan Tissue Explorer[/bold cyan]\n"
            "[dim]Explore glycan profiles, mucins, and lectin-glycan binding.\n"
            "Outputs will auto-save to ./output/[/dim]",
            border_style="cyan", padding=(1, 3),
        )
    )

    MENU = [
        ("1",  "Tissue Exploration (Glycans, Motifs, and Charts)"),
        ("2",  "Mucin Exploration"),
        ("3",  "Mucin Expression Map"),
        ("4",  "Comparing Tissues"),
        ("5",  "Motif Prevalence"),
        ("6",  "Lectin Binding Profile"),
        ("7",  "Lectin & Glycan Binding Strength Comparison"),
        ("8",  "Galectin & Tissue Interaction"),
        ("9",  "Glycan Abundance Ranking in Tissue"),
        ("10", "Binding Score Ranking for a Glycan"),
        ("11", "On-target vs Off-target Glycan Comparison"),
        ("12", "Cross-System Synthesis"),
        ("T",  "Software's Supported Tissues List"),
        ("M",  "Software's Supported Mucins List"),
        ("N",  "Software's Supported Motifs List"),
        ("K",  "Software's Supported Lectins List"),
        ("Q",  "Quit"),
    ]

    while True:
        rprint("")
        rule("Main Menu")
        for key, label in MENU:
            rprint(f"  [cyan]{key:>2}[/cyan]  {label}")
        rprint("")
        
        choice = Prompt.ask("[bold]Choose[/bold]", default="Q").strip().upper()

        if choice == "Q": rprint("[dim]Quitting the Application.[/dim]"); break

        elif choice == "1":
            cmd_list_tissues()
            q = Prompt.ask("Tissue key")
            if q: cmd_cell(q)
        elif choice == "2":
            q = Prompt.ask("Mucin name (e.g. MUC16, CA125)")
            if q: cmd_mucin(q)
        elif choice == "3": cmd_mucin_map()
        elif choice == "4":
            cmd_list_tissues()
            raw = Prompt.ask("Tissue keys separated by spaces")
            if raw: cmd_compare(raw.split())
        elif choice == "5":
            cmd_list_motifs()
            mk = Prompt.ask("Motif key")
            if not mk: continue
            cmd_list_tissues()
            raw = Prompt.ask("Tissue keys separated by spaces")
            if raw: cmd_motif(mk, raw.split())
        elif choice == "6":
            q = Prompt.ask("Lectin (e.g. Gal-3)")
            if q: cmd_binding(q)
        elif choice == "7":
            lec    = Prompt.ask("Lectin (e.g. Gal-3)")
            if not lec: continue
            glycan = Prompt.ask(
                "Glycan — shorthand (e.g. core1, stn) or full IUPAC string"
            )
            if lec and glycan: cmd_binding_pair(lec, glycan)
        elif choice == "8":
            cmd_list_tissues()
            q = Prompt.ask("Tissue key")
            if q: cmd_binding_galectins(q)
        elif choice == "9":
            cmd_list_tissues()
            q = Prompt.ask("Tissue key")
            if q: cmd_tissue_rank(q)
        elif choice == "10":
            g = Prompt.ask("Glycan (shorthand or IUPAC)")
            if g: cmd_binding_rank(g)
        elif choice == "11":
            cmd_list_tissues()
            tissue = Prompt.ask("Tissue key")
            if not tissue: continue
            rprint("\n[dim]Available glycan shorthands:[/dim]")
            for alias, iupac in list(GLYCAN_ALIASES.items()):
                rprint(f"  [cyan]{alias:20s}[/cyan]  [dim]{iupac}[/dim]")
            on_g  = Prompt.ask("On-target glycan (shorthand or IUPAC)")
            off_g = Prompt.ask("Off-target glycan (shorthand or IUPAC)")
            if on_g and off_g:
                cmd_target_compare(tissue, on_g, off_g)
        elif choice == "12":
            cmd_cross_system_synthesis()
        elif choice == "T": cmd_list_tissues()
        elif choice == "M": cmd_list_mucins()
        elif choice == "N": cmd_list_motifs()
        elif choice == "K": cmd_list_lectins()
        else: rprint(f"[red]Unknown option:[/red] {choice!r}")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="python main.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            Glycan Tissue Explorer
            ======================
            All outputs (charts, CSV) auto-save to ./output/

            QUICK START
            -----------
              python main.py                    
              python main.py --cell eye
              python main.py --mucin MUC16
              python main.py --binding Gal-3
              python main.py --target-chart eye
              python main.py --disease-compare dry_eye
              python main.py --cross-system-synthesis
        """),
    )

    listing = parser.add_argument_group("Listing")
    listing.add_argument("--list-systems",  action="store_true")
    listing.add_argument("--list-tissues",  action="store_true")
    listing.add_argument("--list-mucins",   action="store_true")
    listing.add_argument("--list-motifs",   action="store_true")
    listing.add_argument("--list-lectins",  action="store_true")

    analysis = parser.add_argument_group("Analysis")
    analysis.add_argument("--cell",               metavar="TISSUE")
    analysis.add_argument("--mucin",              metavar="NAME")
    analysis.add_argument("--mucin-map",          action="store_true")
    analysis.add_argument("--compare",            nargs="+", metavar="TISSUE")
    analysis.add_argument("--motif",              nargs="+", metavar="KEY_OR_TISSUE")
    analysis.add_argument("--binding",            metavar="LECTIN")
    analysis.add_argument("--binding-pair",       nargs=2, metavar=("LECTIN", "GLYCAN"))
    analysis.add_argument("--binding-galectins",  metavar="TISSUE")
    analysis.add_argument("--tissue-rank",        metavar="TISSUE")
    analysis.add_argument("--binding-rank",       metavar="GLYCAN")
    analysis.add_argument("--target-compare",     nargs=3,
                          metavar=("TISSUE", "ON_GLYCAN", "OFF_GLYCAN"),
                          help="Compare on-target vs off-target glycan in a tissue.  "
                               "E.g. --target-compare eye core1 stn")
    analysis.add_argument("--cross-system-synthesis", action="store_true",
                          help="§ 12 — Cross-system synthesis heatmap")

    opts = parser.add_argument_group("Options")
    opts.add_argument("--top", type=int, default=20, metavar="N")
    return parser


def main(argv=None):
    parser = build_parser()
    args   = parser.parse_args(argv)
    _ensure_output()
    ran = False

    if args.list_systems:    cmd_list_systems();  ran = True
    if args.list_tissues:    cmd_list_tissues();  ran = True
    if args.list_mucins:     cmd_list_mucins();   ran = True
    if args.list_motifs:     cmd_list_motifs();   ran = True
    if args.list_lectins:    cmd_list_lectins();  ran = True
    if args.cell:            cmd_cell(args.cell, top=args.top);  ran = True
    if args.mucin:           cmd_mucin(args.mucin);              ran = True
    if args.mucin_map:       cmd_mucin_map();                    ran = True
    if args.compare:         cmd_compare(args.compare);          ran = True
    if args.motif:
        if len(args.motif) < 2:
            rprint("[red]--motif needs a motif key + at least one tissue.[/red]")
        else:
            cmd_motif(args.motif[0], args.motif[1:])
        ran = True
    if args.binding:         cmd_binding(args.binding);          ran = True
    if args.binding_pair:
        cmd_binding_pair(*args.binding_pair);                    ran = True
    if args.binding_galectins:
        cmd_binding_galectins(args.binding_galectins);           ran = True
    if args.tissue_rank:     cmd_tissue_rank(args.tissue_rank, top=args.top); ran = True
    if args.binding_rank:    cmd_binding_rank(args.binding_rank, top=args.top); ran = True
    if args.target_compare:
        tissue, on_g, off_g = args.target_compare
        cmd_target_compare(tissue, on_g, off_g);                 ran = True
    if args.cross_system_synthesis:
        cmd_cross_system_synthesis();                            ran = True

    if not ran:
        interactive_menu()


if __name__ == "__main__":
    main()