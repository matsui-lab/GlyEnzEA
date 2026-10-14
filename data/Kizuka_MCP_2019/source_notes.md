# Nakano 2019 mouse brain N-glycome (Mgat3 knockout) — source notes

**Citation**: Nakano M, Mishra SK, Tokoro Y, Sato K, Nakajima K, Yamaguchi Y, Taniguchi N, Kizuka Y.
*Bisecting GlcNAc is a general suppressor of terminal modification of N-glycan.* Mol Cell Proteomics 2019;18(10):2044-2057.
doi:10.1074/mcp.RA119.001534. The two workbooks were transcribed by K. Himori from the tables and supplementary data of that paper (quantities and GlyTouCan candidate assignments); they are not files received from the authors.

## Files and what is read from them

| File | Sheet | Used by | Content |
|---|---|---|---|
| `Kizuka_MCP_2019.xlsx` | `quant` | `scripts/mouse/step2_glycome.R` | one row per observed glycan species (`No.` = 1 ... 30 with letter suffixes), peak intensities and relative abundance (%) in Mgat3+/+ and Mgat3−/− brain |
| `20260215_Kizuka_MCP_2019_revise.xlsx` | `anno` | `scripts/mouse/step2_glycome.R` | one row per candidate structure: `No`, monosaccharide composition, `GlyTouCanID_new` (candidate structure; converted to WURCS through `cache/`), `vsWT` (the paper's label: lost / unchanged / increased / control) |

The other sheets (`anno` of the first file, `データベース優先`, `Supple`) are earlier annotation drafts and are not read.

## Species numbering and grouping

- A number with letter suffixes (10a, 10b, 10c, ...) denotes isomers of one composition that the paper lists
  separately. The analysis treats one *observed glycan species* as one abundance measurement, so suffixed rows
  are pooled by number (10a-10c -> one species "10a-c") unless the paper itself keeps them distinct.
- Kept distinct because they are separate measurements in the paper: 9a / 9b, 18a / 18b, 19a / 19bc
  (`nos_keep_separate` in `step2_glycome.R`; quant rows 19b and 19c are joined to the annotation's 19bc).
- Result: 34 observed glycan species, 98 candidate structures; 33 species carry at least one structure and
  enter the ranking (Methods 2.4).

## 10d is kept apart (decision of 2026-10-07)

| species | Mgat3+/+ (%) | Mgat3−/− (%) | structure (GlyTouCan) | paper label |
|---|---|---|---|---|
| 10a | 1.0 | 0.0 | G34428LG / G16146QN (two candidates for 10a-c) | lost |
| 10b | 5.3 | 0.2 | " | lost |
| 10c | 2.2 | 0.1 | " | lost |
| 10d | 0.5 | 2.1 | none | unchanged |

10d has no structure in the paper (no drawing in Fig. 1 or the supplement) and therefore cannot be mapped to any
enzyme. Until 2026-10-07 the code pooled its abundance into "10" together with 10a-10c, so a species that *increases*
in the knockout was added to three that disappear, and the log2 fold change of 10a-c came out as −1.9 instead of
−4.8; the coauthor's Fig. 4A also showed a 10d bar that duplicated the pooled value. Since 2026-10-07 `10d` is in
`nos_keep_separate`: it is its own species, has no structure, and so is excluded from mapping and ranking (the 34th
species that is not ranked).

Effect on the results (all other conclusions unchanged; 16 glycan sets, 33 ranked species, MGAT3 the only set
significant after correction):

| | pooled (v7 numbers) | 10d apart (v8 numbers) |
|---|---|---|
| MGAT3 glycan set NES | −1.66 | −1.60 |
| nominal P | 5.8e-4 | 3.0e-3 |
| Storey q | 0.0062 | 0.033 |
| BH-adjusted P | 0.0093 | 0.049 |
| minimum ES (Fig. 4B) | −0.94 | −0.90 |
| ORA P, curveball P | 7.7e-5, 0.001 | unchanged |

Storey q for MGAT3 stays below 0.05 for every nPermSimple in 10^3-10^5 and every seed (max 0.042); BH does so
in every seed only from nPermSimple = 10^4 (`provenance/log_rerun_20261007_10d.txt`).
