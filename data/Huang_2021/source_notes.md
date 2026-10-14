# Huang 2021 HEK293 N-glycan KO panel — source notes

**File**: `NIHMS1683936-supplement-5.xlsx`
**Citation**: Huang Y-F, Aoki K, Akase S, Ishihara M, Liu Y-S, Yang G, Kizuka Y, Mizumoto S, Tiemeyer M, Gao X-D, Aoki-Kinoshita KF, Fujita M. *Global mapping of glycosylation pathways in human-derived cells.* Dev Cell 2021 Apr 19;56(8):1195-1209.e7. doi:10.1016/j.devcel.2021.02.023 (PMID: 33730547, NIHMSID: NIHMS1683936)

> 2026-08-27 訂正: 以前ここに書かれていた `56:1083-1098.e5` と `PMID: 33765438` は誤り。
> 33765438 は別論文（Molinari M, *ER-phagy responses…*, Dev Cell 56(7):949-966）。
> PubMed esummary で確認して修正した。
**Obtained from**: PubMed Central, author manuscript NIHMS1683936, supplementary file 5 (downloaded 2026-05-19)

## Content used

- **Sheet `Table_S4_1`**, rows 1-61 (top 61 N-glycan rows; rows 62-73 are
  GSL `+Hex2Cer` series and are excluded by glyEnzMap which targets
  N-glycans only).
- **Sheet `Table_S4_2`** is MS2 signature-ion data and is NOT used.

## Per-row layout (Table_S4_1)

| col | content |
|---|---|
| 1 | N-glycan row number (1..73; rows 62..73 are GSL — excluded) |
| 2 | Structural name (e.g. "M3N2", "SA1Gal1N1M3N2F1") |
| 3 | Structure/Composition flag (mostly blank) |
| 4 | **Component ID** — composition-level GlyTouCan (one per row, with ambiguous WURCS having `x` placeholders) |
| 5 | **Structure ID** — structure-level GlyTouCan(s), `;`-separated when multiple isomer candidates exist for the same composition |
| 6..  | Per-condition replicate columns: WT, MGAT1-KO, MGAT2-KO, MGAT4A&4B-KO, MGAT5-KO, T-KO(MANA1&A2&B1), B4GALNT3&4-KO, SLC35C1-KO — each with Exp1 / Exp2 / Average / Error in Rel % |

## Decisions for glyEnzMap ingestion

### Use Structure ID (col 5), not Component ID

The Component ID is composition-level (anomeric/branching info absent, WURCS contains `x` placeholders) and not suitable for structure-level toposeq construction. Structure IDs are biosynthesis-curated structure-level GlyTouCan entries with fully resolved WURCS strings.

### Multi-ID rows → multiple glycomeRes rows, one quant group

When a row lists multiple Structure IDs (e.g. `G56014GC; G84552CD` for M4N2), these are alternative isomer candidates for **one MS peak / one quant observation**. glyEnzMap's existing `group` convention (Kizuka `manual_structures` pattern in `src/scripts/step2_glycome.R` lines 121-213, OR rule documented in `src/scripts/step3_build_pk.R` lines 198-244) handles this case:

- One xlsx row → one `group` (group id = N-glycans row number)
- N Structure IDs → N rows in glycomeRes with the same `group` and same `logFC`, different `toposeq`
- step3 OR rule attributes an enzyme to the group iff any candidate toposeq matches

**No abundance distribution (1/N split) is applied** — the OR collapse at enzyme level is the right semantic.

### Ambiguous LIN handling

A subset of Structure IDs have WURCS with `|` LIN notation indicating truly multi-parent ambiguity (e.g. sialic acid attached to one of several Gal antennae, position uncertain). `build_tree_wurcs` after the 2026-05-19 fix correctly returns `ambiguous=TRUE` for these and they are dropped from the candidate set. If a row's entire Structure ID set becomes ambiguous, that row contributes no toposeq (skipped from glycomeRes for that KO).

Coverage on the 130 unique Structure IDs (after post-fix verification, 2026-05-19):
- 112 / 130 (86.2%) → toposeq OK
- 18 / 130 → ambiguous_skip (truly multi-parent, biologically uncertain)
- Row-level: 44 all-ok + 5 partial-ok + 6 all-fail (out of 55 valid N-glycan rows; 2 rows are header-bleed false positives, 1 row has Structure ID = "N/A" in the source)

### Quant aggregation

For each row × KO: `logFC = log2((KO_Average + eps) / (WT_Average + eps))` with `eps = 1e-3` (matches Kizuka convention in `src/scripts/step2_glycome.R` line 16). Per-glycan p-value is not computed (n=2 replicates / condition — Kizuka convention is to feed logFC ranks into GSEA and rely on GSEA enrichment p).

## Out-of-scope

- O-glycan and glycoproteomic-site supplementary tables in other supps
  (NIHMS1683936-supplement-2, 3, 4, 6, 7, 8, 9) — N-glycan panel only.
- Quantitative comparisons across KOs other than vs WT (all comparisons
  are KO vs WT, one at a time).
- Use of MS2 signature ions from Table_S4_2.

## Files in this directory

```
NIHMS1683936-supplement-5.xlsx   (248K, the original supp file)
source_notes.md                  (this file)
```
