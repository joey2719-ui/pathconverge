# pathconverge

An R package for testing **panel-level pathway convergence** between a nominated biomarker panel and an independent molecular signature.

Biomarker panels are validated statistically (discrimination, calibration, decision-curve utility), but the biological coherence of a panel with independent molecular layers is never quantified. `pathconverge` asks the panel-level question: *does this panel's pathway annotation converge with an independent molecular signature more often than expected by chance?*

## How it works

- **Signature-permutation null** that preserves panel size, the measured compound universe, and — optionally — feature direction.
- **Specificity-stratified convergence curve** T(k) across signature depths k.
- **Random-panel benchmark** that calibrates the observed convergence against size-matched panels drawn from the same pathway universe, with a primary specification conditioned on measurable links and pathway size, plus uniform and unconditioned sensitivity variants.
- **Dual mode**: metabolite mode (KEGG compound membership) and gene mode (KEGG gene membership); binary co-membership and direction-aware scoring.

## Installation

```r
# from a local clone
install.packages("pathconverge", repos = NULL, type = "source")
# or after GitHub release
# devtools::install_github("<username>/pathconverge")
```

Requires R >= 4.1; imports `data.table` (ggplot2/knitr/rmarkdown optional for plots and the vignette).

## Quick start

```r
library(pathconverge)

# membership: pathway x feature (KEGG compound or gene) 0/1 matrix
# measured:   features quantified in the signature layer
# links:      panel protein -> predicted pathway data.table
# signature:  ranked feature list from the independent molecular layer

res <- permutation_null(links, membership, measured, signature,
                        n_perm = 10000, seed = 20261028)
res$p_value

curve <- convergence_curve(links, membership, measured,
                           ks = c(5, 8, 10, 15, 20, 25, 31), seed = 20261028)
plot_convergence_curve(curve)

bench <- random_panel_benchmark(links, membership, measured, signature,
                                n_panels = 100, seed = 20261028)
```

See `vignettes/nsclc_case_study.Rmd` for the full non-small-cell lung cancer (NSCLC) plasma proteomics–metabolomics case study.

## Analysis pipeline (`analysis/`)

Numbered scripts reproduce the manuscript analyses in run order. Case-study inputs are cached under a local `cache/` directory so repeated runs are offline.

| Script | Purpose |
|---|---|
| `01_kegg_membership.R` | Build the KEGG compound → pathway membership structure (cached, offline reuse) |
| `02_null_model.R` | Seed implementation: permutation null + specificity curve + random-panel benchmark for the NSCLC four-protein case study |
| `03_sensitivity_ora_thresholds.R` | Support-set sensitivity to the ORA significance threshold θ (hypergeometric, KEGG background N = 4702) |
| `04_channel_inference.R` | TMT plex channel inference: are channels 128–131 four stage pools of one pool design? |
| `05_package_case_study.R` | Full case study through the installed package; verifies agreement with seed scripts (01–03) |
| `06_benchmark_universe_check.R` | Benchmark null-specification variants (primary / uniform full-map / uniform metabolic-only / unconditioned) |
| `07_simulation_study.R` | Type-I error, p-value uniformity, and power across 16 configurations anchored on the empirical KEGG compound–pathway bipartite structure |
| `08_hmdb_mapping_expansion.R` | Multi-route HMDB → KEGG mapping and mapping-expansion sensitivity |
| `09_mapping_diagnostics.R` | Mapping diagnostics on the cleaned metabolomics feature table |
| `10_direction_weakening_diag.R` | Cross-test decomposing the direction-preserving p shift (universe vs direction-pool effects) |
| `11_benchmark_competitors.R` | ORA / Fisher / naive competitors on the same simulated configurations, exact RNG-stream pairing |
| `12_benchmark_gsea_ranked.R` | Ranked (fgsea-style) competitor on the same configurations, exact RNG-stream pairing |
| `13_simulation_figure.R` | Two-panel simulation figure (null calibration + power) |
| `14_external_case_gse75037.R` | External case A (gene mode): GSE75037, 83 paired LUAD / non-malignant tissues (GPL6884 annotation) |
| `15_external_case_st000368.R` | External case B (compound mode): ST000368 plasma GC–TOF metabolomics (Metabolomics Workbench PR000293) |
| `16_external_cases_figure.R` | Two-panel external-case figure (GSE75037 gene mode; ST000368 compound mode with the four-protein replication test and six-gene comparison) |
| `17_manuscript_fig3.R` | Combined manuscript Figure 4 (convergence curve + null distribution), with a regression guard on the canonical p-value |
| `18_supplementary_tables.R` | Assemble supplementary tables S1–S4 |
| `20_fourprotein_on_st000368.R` | Replication test: the four-protein discovery panel tested against the independent ST000368 compound-mode cohort |
| `21_link_count_diagnostics.R` | Link-accounting diagnostics: measurable and supported link counts for the external-cohort analyses |
| `22_benchmark_calibration.R` | Type-I error calibration of the random-panel benchmark under the L = 0 simulation configurations |

All scripts use a fixed random seed (`20261028`) for reproducibility.

## Data

- **Case study**: in-house NSCLC plasma proteomics (TMT, 4-protein panel nomination) and metabolomics (KEGG-mapped feature table).
- **External cohorts**: [GSE75037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75037) (Gene Expression Omnibus) and [ST000368](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000368) (Metabolomics Workbench, study PR000293).

## Citation

If you use `pathconverge`, please cite: Lin, C.-H. et al. *pathconverge: a panel-level permutation null for quantifying pathway convergence between biomarker panels and metabolic signatures.* (manuscript in submission; DOI to be minted on Zenodo).

## License

MIT — see [LICENSE](LICENSE).
