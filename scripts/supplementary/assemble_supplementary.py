#!/usr/bin/env python3
"""Assemble paper/supplementary/GlyEnzEA_v8_supplementary.tex from the built tables and figures.
Run from the repository root after build_tables.R and figS1_column_examples.R (+ Chrome print of figS1.html).
Text of the Supplementary Methods lives in this file."""
import csv, pathlib, re
B = pathlib.Path("paper/supplementary/build"); OUT = pathlib.Path("paper/supplementary/GlyEnzEA_v8_supplementary.tex")
def esc(s):
    return s.replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace("#", "\\#")
def brk(s):
    return esc(s).replace("/", "/\\allowbreak{}")

def num(s):
    s = s.strip()
    m = re.fullmatch(r"(-?\d(?:\.\d+)?)e(-?\d+)", s)
    if m: return f"${m.group(1)} \\times 10^{{{int(m.group(2))}}}$"
    return esc(s)
def longtable(path, colspec, caption, label, font="\\scriptsize", landscape=False, numeric_cols=()):
    rows = list(csv.reader(open(path, encoding="utf-8"))); head, body = rows[0], rows[1:]
    hdr = " & ".join(f"\\textbf{{{esc(h)}}}" for h in head) + " \\\\"
    lines = [" & ".join(num(c) if i in numeric_cols else brk(c) for i, c in enumerate(r)) + " \\\\" for r in body]
    t = (f"{font}\n\\begin{{longtable}}{{{colspec}}}\n\\caption{{{caption}}}\\label{{{label}}}\\\\\n\\toprule\n{hdr}\n\\midrule\n\\endfirsthead\n"
         f"\\caption[]{{(continued)}}\\\\\n\\toprule\n{hdr}\n\\midrule\n\\endhead\n\\midrule\\multicolumn{{{len(head)}}}{{r}}{{\\textit{{continued on next page}}}}\\\\\n\\endfoot\n\\bottomrule\n\\endlastfoot\n"
         + "\n".join(lines) + "\n\\end{longtable}\n\\normalsize\n")
    return ("\\newgeometry{textwidth=261mm,textheight=170mm,left=18mm,top=20mm,bottom=20mm}\\special{papersize=297mm,210mm}\n\\setlength{\\tabcolsep}{3pt}\\setlength{\\LTcapwidth}{\\textwidth}\n" + t + "\\restoregeometry\\special{papersize=210mm,297mm}\\setlength{\\LTcapwidth}{170mm}\\setlength{\\tabcolsep}{6pt}\n") if landscape else t

methods = r"""
\section*{Supplementary Methods}\label{sm}

\subsection*{S1. Mouse glycome input and measured features}
The Mgat3-knockout mouse brain N-glycome (Nakano et al., 2019) was used as received from the Kizuka laboratory: one row per
measured features, numbered 1--30 with letter suffixes for isomers of one composition, with relative abundances (\%)
in Mgat3$^{+/+}$ and Mgat3$^{-/-}$ brain and, for most features, one or more candidate structures (GlyTouCan identifiers).
Suffixed rows of one composition were pooled into one measured feature unless the paper itself reports them as
separate measurements (9a/9b, 18a/18b, 19a/19bc). Species 10d, which has no structure in the paper and increases in the
knockout while 10a--10c disappear, is kept apart from 10a--c; having no structure, it is not mapped and not ranked. This
gives 34 measured features with 98 candidate structures, of which 33 features enter the ranking. The relative
abundances are normalized to the total of oligomannose-type N-glycans and are therefore not a closed composition; scores
are $x_i = \log_2[(\bar y_{i,1}+\epsilon)/(\bar y_{i,0}+\epsilon)]$ with $\epsilon = 10^{-3}$ percentage points, one score
per measured feature regardless of the number of candidates.

\subsection*{S2. Filtering order and missing values}
Structural matching (Methods 2.2) first produces a binary matrix of enzymes by candidate topologies (93 enzymes $\times$ 97
topologies in the mouse dataset). Filtering is applied in the following order. (i) At the candidate level, enzymes that
match every candidate topology are removed, then enzymes that match none. (ii) Candidate columns are combined into
measured features by OR, i.e.\ a feature is in an enzyme's glycan set if any of its candidates matches. (iii) At the
feature level, feature columns to which no enzyme maps are removed, then enzymes that map to every remaining feature,
then enzymes that map to none. Both levels are filtered because a rule can match all candidates without matching all
feature, and vice versa. In the mouse dataset this leaves 48 enzymes and 33 features; identical rows are then collapsed
into the 16 distinct glycan sets of the hypothesis family. Missing values arise only from features without a structure
(10d), which are dropped before matching; scores are never missing because of the pseudocount. The ranked universe is the
set of features with at least one structure and a finite score (33), and features are not merged when their enzyme-mapping
patterns coincide.

\subsection*{S3. Reaction rules and matching criteria}
Human and mouse N-glycan reactions were taken from GlycoEnzDB and GlyGen and represented as enzyme--reactant--product
triplets. Carrier prefixes and notation-specific decorations were normalized, reactions with unresolved sulfate or
phosphate placement were excluded, and canonical topologies were deduplicated, giving 2,316 triplets for 93 enzymes.
Whether a rule $r = (e, s, p)$ matches an observed structure $g$ is decided by subgraph containment ($\preceq$, VF2 in
igraph) of the canonical topologies, with the criterion depending on the source and class of the rule:

\begin{center}\small
\begin{tabular}{lll}
\toprule
Source & Enzyme class & Criterion \\
\midrule
GlycoEnzDB & Glycosyltransferase & $\phi(s)\preceq\phi(g)$ and $\phi(p)\preceq\phi(g)$ \\
GlycoEnzDB & Glycosidase & $\phi(p)\preceq\phi(g)$ and not $\phi(s)\preceq\phi(g)$ \\
GlycoEnzDB & (empty) or N\_linked & as Glycosyltransferase (joint reactant--product criterion) \\
GlyGen & all & as Glycosyltransferase (joint reactant--product criterion) \\
\bottomrule
\end{tabular}
\end{center}

Reactant and product containment are evaluated independently, without a shared reaction-centre alignment; neither
root-to-root alignment nor compartment is imposed. Matches from several rules of one enzyme are combined by OR.

\subsection*{S4. Enrichment, correction and null checks}
ORA uses one-sided Fisher tests on the source-annotated increased and lost lists against the union of valid annotated
features and set members. Preranked GSEA uses fgsea (adaptive multilevel estimation, nPermSimple = 50,000, minSize = 1,
maxSize = 500 for the mouse dataset and no upper limit for HEK293, eps = 0 and seed 42 for HEK293). The primary correction
is the Storey q-value (qvalue, bootstrap estimation of $\pi_0$) within each hypothesis family (16 sets in the mouse
dataset, 22 per HEK293 perturbation); when $\pi_0$ cannot be estimated the Benjamini--Hochberg adjustment is used
instead, and Tables S1 and S2 state which was used. For the mouse dataset Storey estimation succeeded; for HEK293 it
succeeded in two of the seven perturbations (MGAT2 and B4GALNT3/4). Three null checks accompany the mouse result:
1,000 degree-preserving curveball randomizations of the enzyme--feature matrix (empirical resolution 0.001), reassignment
of the 16 set labels (resolution 1/16), and size-matched random sets. Sensitivity to minSize (1--3) and to nPermSimple
($10^3$--$10^5$, five seeds each) is reported in Table S1 and Figure S2.

\subsection*{S5. HEK293 perturbation panel}
From the N-glycan section of Table S4 of Huang et al.\ (2021), 41 rows were taken and one row without a structure
identifier was excluded, leaving 40 measured features. The 107 GlyTouCan identifiers assigned to these rows were
converted to WURCS (107/107) and canonicalized as in Methods 2.2. Two sources of multiple candidates were kept rather than
resolved: 24 rows carry several identifiers (up to 13), and 30 of the 107 WURCS strings list alternative linkages, of
which 18 leave the parent residue ambiguous and were expanded into all concrete structures before canonicalization
(e.g.\ 729 combinations reduce to 304 topologies, 256 to 13, 8 to 2). This gives 823 feature--candidate pairs (822
distinct topologies). Scores are $\log_2$ ratios of the reported condition means with the same $\epsilon$ as above.
One mapping (60 enzymes, 22 distinct glycan sets, set sizes 1--29) serves all seven perturbations; only the ranking
differs. Target sets and expected directions were listed before the analysis from the deleted reactions; SLC35C1, a
GDP-fucose transporter without an enzyme rule, is read through the FUT8 glycan set. No correction is pooled across
perturbations.
"""

capS1 = ("Candidate structures in each display column of Figure 3. One row per column of Figure 3 (G labels), in the "
         "order of Figure 3 from left to right; each row shows the candidate structures whose topology has that "
         "enzyme-mapping pattern, drawn in SNFG without linkage labels, with the measured features (numbering of "
         "Nakano et al., 2019) to which each candidate belongs. Several candidate structures can share one topology "
         "because linkage positions are not retained, so a row can show more drawings than topologies. Candidates of "
         "features 6, 26ab and 30 (asterisks) were entered manually without GlyTouCan identifiers and are drawn from "
         "their topology.")
capS2 = ("Null checks for the mouse enrichment results. (A) fgsea score-permutation $P$ against the empirical $P$ "
         "from 1,000 degree-preserving curveball randomizations of the enzyme--feature matrix, both $-\\log_{10}$; "
         "point size is set size, colour the direction of NES; dashed line $y=x$, dotted lines $P = 0.05$. The MGAT3 "
         "glycan set is the most extreme set under both nulls and reaches the empirical resolution of the curveball "
         "null ($P = 0.001$). (B) fgsea $P$ (mean $\\pm$ SD over five seeds) for nPermSimple from $10^3$ to $10^5$; the "
         "dashed line marks $P = 0.05$. Point estimates are stable; for the MGAT3 glycan set the Storey q-value is "
         "below 0.05 in every run.")

import subprocess
def figS1_pages(cap):
    n = int([l for l in subprocess.run(["pdfinfo", str(B / "figS1.pdf")], capture_output=True, text=True).stdout.splitlines() if l.startswith("Pages:")][0].split()[1])
    out = ["\\begin{figure}[p]\\centering\\includegraphics[page=1,width=\\textwidth]{build/figS1.pdf}\\caption{" + cap + f" (page 1 of {n})" + "}\\label{figS1}\\end{figure}"]
    for k in range(2, n + 1):
        out.append("\\clearpage\\begin{figure}[p]\\centering\\includegraphics[page=%d,width=\\textwidth]{build/figS1.pdf}\\caption*{\\textbf{Figure S1} (continued, page %d of %d).}\\end{figure}" % (k, k, n))
    return "\n".join(out)

tex = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=20mm]{geometry}
\usepackage{amsmath,amssymb,booktabs,longtable,graphicx,caption,array}
\captionsetup{font=small,labelfont=bf,justification=raggedright,singlelinecheck=false}
\setlength{\parskip}{4pt}\setlength{\parindent}{0pt}\setlength{\LTcapwidth}{170mm}
\renewcommand{\thefigure}{S\arabic{figure}}\renewcommand{\thetable}{S\arabic{table}}
\begin{document}
\begin{center}{\large\textbf{Supplementary Information}}\\[4pt]
GlyEnzEA: reaction-informed enrichment of enzyme-centered glycan sets from glycomics data\\
Koichi Himori, Bingyuan Zhang, Yusuke Matsui\end{center}
\vspace{4pt}
\begin{tabular}{@{}p{150mm}l@{}}
Supplementary Methods S1--S5 & page \pageref{sm} \\
Table S1. Glycan sets of the Mgat3-knockout mouse glycome: composition, GSEA, correction, null checks \\quad & page \pageref{tabS1} \\
Table S2. GSEA of all 22 glycan sets in the seven HEK293 perturbations & page \pageref{tabS2} \\
Figure S1. Candidate structures in each display column of Figure 3 & page \pageref{figS1} \\
Figure S2. Null checks for the mouse enrichment results & page \pageref{figS2} \\
\end{tabular}
""" + methods + r"""
\clearpage
""" + longtable(B / "TableS1_mouse_sets.csv", "l>{\\raggedright\\arraybackslash}p{3.6cm}lrrrrrlrrrr",
      "Enzyme-centered glycan sets of the Mgat3-knockout mouse glycome: composition, preranked GSEA, correction, "
      "degree-preserving (curveball) null and minSize sensitivity. NES and nominal $P$ from fgsea; Storey q is the "
      "primary correction and BH the Benjamini--Hochberg adjustment within the family of 16 sets; curveball values are "
      "from 1,000 degree-preserving randomizations; the last two columns are the fgsea $P$ with minSize 2 and 3 (sets "
      "below the threshold are excluded, shown as NA).", "tabS1", font="\\scriptsize", landscape=True, numeric_cols=(5, 6, 7, 10, 11, 12)) + r"""
\clearpage\setlength{\tabcolsep}{3.5pt}
""" + longtable(B / "TableS2_hek293_all.csv", ">{\\raggedright\\arraybackslash}p{2.4cm}>{\\raggedright\\arraybackslash}p{3.2cm}crrrrrl",
      "Preranked GSEA of all 22 enzyme-centered glycan sets in each of the seven HEK293 perturbations (Huang et al., "
      "2021). Target set marks the set expected to change from the deleted reactions. q is the within-perturbation "
      "value used in the paper: the Storey q-value, or the BH adjustment (BH*) when $\\pi_0$ could not be estimated, as "
      "stated in the last column; BH $P$ is the BH-adjusted $P$ for every set.", "tabS2", font="\\scriptsize", numeric_cols=(5, 6, 7)) + r"""
\clearpage
""" + figS1_pages(capS1) + r"""
\clearpage
\begin{figure}[p]\centering\includegraphics[width=\textwidth]{../../output/FigureS_stat_robustness.pdf}
\caption{""" + capS2 + r"""}\label{figS2}\end{figure}
\end{document}
"""
OUT.write_text(tex, encoding="utf-8"); print("wrote", OUT)
