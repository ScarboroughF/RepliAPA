# RepliAPA

**Replicate-aware compositional inference of alternative polyadenylation in
single-cell perturbation screens**

RepliAPA is an R package for target-level inference of perturbation-associated
alternative polyadenylation (APA). It aggregates cell measurements to
guide-by-biological-replicate pseudobulks, retains whole guides as intervention
units, and jointly evaluates all poly(A) sites (PASs) assigned to a gene.

RepliAPA starts from an existing PAS-by-cell count matrix. It does **not** read
FASTQ files, detect PASs, perform guide assignment, or treat cells as
independent biological replicates.

> **Release status:** version 0.99.0 is a pre-release candidate in the public
> repository https://github.com/ScarboroughF/RepliAPA. The archived source package
> previously passed `R CMD check --no-manual` with 0 errors, 0 warnings and 0 notes.
> The versioned release is `v0.99.0`; an archive DOI remains pending.

## Why RepliAPA?

Single-cell perturbation screens may contain thousands of cells but only a few
independent guides and biological replicates. Cells assigned to the same guide
share an intervention and therefore do not constitute independent perturbation
replicates. RepliAPA uses the following inferential hierarchy:

```text
cell measurement
    -> guide x biological-replicate pseudobulk
    -> guide-level intervention evidence
    -> perturbation-target effect
    -> target x gene APA inference
```

For target `t` and gene `g`, the primary estimand is the target-associated
change in the complete PAS composition vector relative to non-targeting
controls (NTCs). RepliAPA reports the gene-global test, the PAS-effect vector,
guide-level evidence, guide heterogeneity, coordinate-based length shift when
positions are available, filtering and explicit abstention.

## Position in an analysis workflow

```text
FASTQ / aligned reads
    -> cell and guide QC
    -> PAS detection and quantification
       (for example PASTA, SCAPE, scTail, Sierra or another caller)
    -> fixed PAS-by-cell count matrix
    -> RepliAPA pseudobulk and inference
    -> statistical results, diagnostics and figures
```

PAS detection and caller comparison are upstream tasks. RepliAPA does not alter
the supplied site set to improve downstream results.

## Main capabilities

- matrix, sparse `Matrix`, Seurat, `SingleCellExperiment` and
  `SummarizedExperiment` input;
- count-conserving cell-to-guide-by-replicate pseudobulk aggregation;
- guide-aware finite-sample permutation-calibrated all-PAS inference;
- target-wise multi-gene fitting with unsupported-target accounting;
- PAS-level target effects and guide-by-replicate effects;
- guide heterogeneity and same-target guide reproducibility;
- whole-guide NTC pseudo-target calibration;
- independently fitted biological-replicate comparison;
- matched all-PAS versus proximal/distal comparison;
- APA and total gene-expression effects reported separately;
- cross-context and post-freeze orthogonal comparison with adequacy guards;
- DRIMSeq input export and matched/native tested-universe accounting;
- deterministic result bundles with SHA256 verification;
- publication-oriented visualizations returned as editable `ggplot2` objects;
- explicit tested, filtered and abstained hypothesis accounting.

## Requirements

- R >= 4.2.0
- imported packages: `digest`, `Matrix`, `methods`, `stats`, `utils`
- `ggplot2` >= 3.4.0 for visualization
- optional container support through `SeuratObject`,
  `SingleCellExperiment` or `SummarizedExperiment`

RepliAPA is implemented in R and does not require a compiled backend.

## Installation

### From GitHub

The public repository is available at https://github.com/ScarboroughF/RepliAPA.

```r
install.packages("remotes")
remotes::install_github(
  "ScarboroughF/RepliAPA",
  build_vignettes = TRUE
)
```

For manuscript analyses, install an immutable tagged release rather than the
moving default branch:

```r
remotes::install_github(
  "ScarboroughF/RepliAPA@v0.99.0",
  build_vignettes = TRUE
)
```

Release notes and the installable R source archive are provided at
https://github.com/ScarboroughF/RepliAPA/releases/tag/v0.99.0.

### From the checked source archive

```r
install.packages(
  "RepliAPA_0.99.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

Confirm the installed version and citation metadata:

```r
library(RepliAPA)
packageVersion("RepliAPA")
citation("RepliAPA")
```

## Input contract

RepliAPA requires three aligned objects.

### 1. PAS count matrix

- rows: unique PAS identifiers;
- columns: cells before aggregation, or guide-by-replicate sample units after
  aggregation;
- values: finite, non-negative integer counts;
- row and column names: required and unique.

### 2. Cell or sample metadata

Row names must match count-matrix column names. For cell-level aggregation the
metadata must contain:

| Field | Meaning |
|---|---|
| `cell_line` | cell line or biological context |
| `replicate` | independent biological replicate |
| `guide` | guide identity |
| `target` | perturbation target; NTCs may use `NT` or another declared label |
| `ntc` | logical non-targeting-control indicator |

Guide-assignment confidence, multiplet removal and cell QC should be resolved
upstream. RepliAPA does not infer missing guide or replicate labels.

### 3. PAS metadata

Row names must match count-matrix row names. The metadata must contain:

| Field | Requirement | Meaning |
|---|---|---|
| `gene` | required | unique PAS-to-gene assignment |
| `position` | optional | genomic PAS coordinate |
| `strand` | optional | `+` or `-` |

`position` and `strand` are required for coordinate-direction and expected
length-shift summaries, but not for the all-PAS global test.

## Quick start with deterministic test data

The bundled example contains two biological replicates, three targeting
guides, four NTC guides and three genes. It demonstrates the interface only and
must not be used to support performance claims.

```r
library(RepliAPA)

x <- repliapa_example()

input <- repliapa_read_pas_counts(
  x$counts,
  cell_data = x$sample_data,
  pas_data = x$pas_data
)

fit <- repliapa_fit(
  input$counts,
  input$cell_data,
  input$pas_data,
  target = "T1"
)

fit
fit$gene_table
fit$pas_effects
fit$guide_effects
fit$heterogeneity
```

## Cell-level workflow

If the columns are individual cells, aggregate them before fitting:

```r
input <- repliapa_read_pas_counts(
  pas_by_cell_counts,
  cell_data = cell_metadata,
  pas_data = pas_metadata
)

pb <- repliapa_pseudobulk(
  input$counts,
  input$cell_data,
  cell_line = "cell_line",
  replicate = "replicate",
  guide = "guide",
  target = "target",
  ntc = "ntc"
)

stopifnot(sum(pb$counts) == sum(input$counts))

fit <- repliapa_fit(
  pb$counts,
  pb$sample_data,
  input$pas_data,
  target = "NUDT21"
)
```

The columns of `pb$counts`, not the original cells, are the sample units passed
to the inferential engine.

## Seurat and Bioconductor containers

Container input requires an explicit PAS assay. This fail-closed behavior
prevents the default RNA assay from being analyzed accidentally.

```r
seurat_input <- repliapa_read_pas_counts(
  seurat_object,
  assay = "polyA",
  layer = "counts"
)

sce_input <- repliapa_read_pas_counts(
  sce_object,
  assay = "polyA"
)

se_input <- repliapa_read_pas_counts(
  summarized_experiment,
  assay = "polyA"
)
```

## Pre-fit adequacy audit

Filtering and support thresholds should be specified using development data and
frozen before validation data are inspected.

```r
audit <- repliapa_audit_input(
  pb$counts,
  pb$sample_data,
  input$pas_data,
  development_replicates = "R1",
  min_cells = 25L,
  min_polya_count = 20L,
  min_target_guides = 2L
)

audit$summary
audit$contract_checks
audit$sample_units
audit$target_support
audit$pas_adequacy
audit$gene_adequacy
```

Thresholds above are illustrative. They must be justified and reported for
each study rather than tuned on frozen validation outcomes.

## Primary model and important parameters

```r
fit <- repliapa_fit(
  pb$counts,
  pb$sample_data,
  input$pas_data,
  target = "NUDT21",
  gene = "gene",
  replicate = "replicate",
  guide = "guide",
  target_col = "target",
  ntc = "ntc",
  replicates = NULL,
  pseudocount = 0.5,
  min_gene_total = 20,
  max_combinations = 1000000
)
```

| Parameter | Purpose |
|---|---|
| `target` | perturbation target tested against NTC guides |
| `replicates` | predeclared biological replicates to include |
| `pseudocount` | positive value used before composition transformation |
| `min_gene_total` | minimum total PAS count in the tested contrast |
| `max_combinations` | maximum admissible guide allocations enumerated |

At least two complete target guides and two complete NTC guides are required.
When the admissible allocation count exceeds `max_combinations`, RepliAPA
abstains explicitly; it does not silently split guides or sample cell labels.
Complete argument documentation is available from R:

```r
?repliapa_fit
?repliapa_audit_input
?repliapa_ntc_calibration
```

## Result objects

`fit$gene_table` contains one row per gene:

| Column | Interpretation |
|---|---|
| `p_value` | finite guide-allocation gene-global APA P value |
| `q_value` | BH-adjusted value within the fitted target |
| `effect_norm` | magnitude of the complete PAS composition change |
| `coordinate_shift` | position-aware composition shift, when available |
| `permutation_n` | number of admissible guide allocations |
| `tested` | gene entered the statistical test |
| `filtered` | gene was excluded by a declared adequacy rule |
| `abstained` | model declined to test despite otherwise available input |
| `reason` | machine-readable filtering or abstention reason |

Additional components are:

- `fit$pas_effects`: target-level effect for each PAS, including replicate
  effects;
- `fit$guide_effects`: each guide-by-replicate PAS effect relative to NTCs;
- `fit$heterogeneity`: within-target guide similarity and dispersion.

A gene-global P value does not imply that every PAS component is individually
significant. Statistical evidence, effect magnitude, guide agreement,
biological-replicate support and tested-universe accounting should be reported
together.

## Null calibration

NTC pseudo-target contrasts relabel whole guides while retaining their
replicate blocks. Cells from one guide are never divided across pseudo-groups.

```r
null <- repliapa_ntc_calibration(
  pb$counts,
  pb$sample_data,
  input$pas_data,
  guide_sizes = 2L,
  n_contrasts = 500L,
  seed = 20260819L
)

null$summary
repliapa_plot_calibration(null)
```

Precomputed contrasts and seeds should be archived for manuscript analyses.
The word “exact” should not be used unless the experimental randomization
distribution and its assumptions have been established independently.

## Biological-replicate validation

Development and validation replicates must be fitted separately. Random cell
splits are not biological-replicate validation.

```r
development <- repliapa_fit(
  pb$counts, pb$sample_data, input$pas_data,
  target = "NUDT21", replicates = "R1"
)

validation <- repliapa_fit(
  pb$counts, pb$sample_data, input$pas_data,
  target = "NUDT21", replicates = "R2"
)

replication <- repliapa_reproducibility(development, validation)
attr(replication, "summary")
repliapa_plot_reproducibility(replication)
```

## All-PAS versus proximal/distal representation

The same guide units and testing engine are used so that only the PAS
representation changes.

```r
pd_fit <- repliapa_fit_proximal_distal(
  pb$counts,
  pb$sample_data,
  input$pas_data,
  target = "NUDT21"
)

representation <- repliapa_compare_allpas_pd(fit, pd_fit)
attr(representation, "accounting")
repliapa_plot_allpas(representation)
```

All-PAS-only discoveries still require independent replicate and guide support.

## APA and total gene expression

PAS composition and total gene-expression effects are separate estimands.

```r
ge <- repliapa_gene_expression_effect(
  gene_expression_counts,
  pb$sample_data,
  target = "NUDT21"
)

classes <- repliapa_classify_effects(fit, ge)
repliapa_plot_apa_ge(classes)
```

RepliAPA does not impute PAS usage from total gene expression.

## Baseline interoperability

Export the unfiltered guide-by-replicate contrast for DRIMSeq:

```r
drim_input <- repliapa_drimseq_input(
  pb$counts,
  pb$sample_data,
  input$pas_data,
  target = "NUDT21"
)

drim_input$counts
drim_input$samples
drim_input$accounting
```

Comparator results should be evaluated on both the method-native tested set and
a common matched tested-gene set. Tested-gene, filtering, abstention, runtime
and memory accounting must not be omitted.

## Publication-quality visualization

```r
p <- repliapa_plot_pas_effect(
  fit,
  gene = "GENE_A",
  pas_data = input$pas_data
)

p <- p + ggplot2::labs(
  title = "Target-associated multi-PAS composition shift"
)

ggplot2::ggsave(
  "pas_effect.pdf",
  p,
  width = 8,
  height = 5,
  device = grDevices::cairo_pdf
)

ggplot2::ggsave(
  "pas_effect.png",
  p,
  width = 8,
  height = 5,
  dpi = 360
)
```

Available result-native plots include:

- `repliapa_plot_pas_effect()`;
- `repliapa_plot_guide_evidence()`;
- `repliapa_plot_calibration()`;
- `repliapa_plot_reproducibility()`;
- `repliapa_plot_allpas()`;
- `repliapa_plot_context()`;
- `repliapa_plot_apa_ge()`.

Each function returns an editable `ggplot` object. Plotting does not refit,
filter or reinterpret the statistical result.

## Reproducible result bundles

```r
bundle <- repliapa_result_bundle(fit, input_audit = audit)
repliapa_verify_contract(bundle = bundle)

repliapa_write_bundle(bundle, "repliapa_results")
repliapa_verify_bundle("repliapa_results")

sessionInfo()
```

The bundle records tested, filtered and abstained hypotheses and includes
SHA256 checksums. Runtime and memory remain `NOT_INSTRUMENTED` unless measured
externally with an explicitly recorded scope.

## Documentation and test data

- `vignettes/end-to-end-repliapa.Rmd`: executable workflow from standardized
  input through verified results;
- `vignettes/visualizing-repliapa-results.Rmd`: editable publication figures;
- `repliapa_example()`: deterministic interface test data;
- `tests/testthat/`: input, inference, governance, validation, visualization and
  API-contract tests;
- `NEWS.md`: versioned changes;
- `DEVELOPMENT_STATUS.md`: candidate integrity and compatibility history.

After installation:

```r
browseVignettes("RepliAPA")
help(package = "RepliAPA")
```

## Validation scope and limitations

RepliAPA has been evaluated using simulation, real NTC semi-simulation,
biological-replicate validation, same-target guide comparisons, a second cell
context from the same parent study and post-freeze MPRA direction support.
These evaluations do not establish universal superiority over sample-aware
alternatives.

Important limitations are:

- finite guide-allocation P values may be discrete when few guides are
  available;
- at least two complete target guides and two complete NTC guides are required;
- biological-replicate support depends on independently generated replicates;
- the package does not estimate a population distribution of heterogeneity
  from only two replicates;
- PAS detection quality and PAS-to-gene mapping remain upstream dependencies;
- cross-cell-line evidence is not patient-level or multicentre validation;
- unsupported targets and genes should remain visible in reporting.

The strongest intended interpretation is a calibration–replicate-supported
discovery trade-off with retained guide-level and internal-PAS information, not
universal dominance over DRIMSeq or other methods.

## Data availability

The method-development analyses use public CPA-Perturb-seq data from GEO
SuperSeries `GSE269600`, including scRNA-seq SubSeries `GSE269596` and MPRA
SubSeries `GSE269595`. RepliAPA does not redistribute those source datasets.
Users should obtain them from the original repositories and comply with the
associated terms and metadata.

The exact manuscript software release, test data, analysis scripts, frozen
configuration and figure source data must be deposited in a DOI-issuing archive
before submission:

```text
Source repository: https://github.com/ScarboroughF/RepliAPA
Versioned release: https://github.com/ScarboroughF/RepliAPA/releases/tag/v0.99.0
Archive DOI:       PENDING
```

## Citation

Until the method paper and archive DOI are available, cite the software using:

```r
citation("RepliAPA")
```

Current software citation:

> Xu W. RepliAPA: Replicate-Aware Inference for Alternative
> Polyadenylation. R package version 0.99.0, 2026.

The repository also contains `CITATION.cff` for automated citation discovery.
The final paper DOI and archived software DOI must be added after they exist;
they must not be anticipated or fabricated.

## Reporting checklist

When reporting a RepliAPA analysis, include:

1. upstream PAS caller and fixed site-set definition;
2. cell, guide, target and biological-replicate counts;
3. targeting- and NTC-guide adequacy;
4. pseudobulk definition and count-conservation check;
5. predeclared filtering and abstention rules;
6. complete tested-gene universe;
7. gene-global P values, multiple-testing procedure and PAS effects;
8. guide and biological-replicate support;
9. APA and gene-expression effects separately;
10. matched and method-native baseline comparisons;
11. seeds, software versions, runtime and memory scope;
12. source-code release, immutable archive and data accessions.

## Support

For reproducible bug reports, include `sessionInfo()`, the RepliAPA version, a
minimal example, expected behavior and complete error output.

- GitHub issues: available after repository publication;
- maintainer: Wei Xu (徐炜), `xw222300@163.com`.

Please do not send controlled, identifiable or unpublished patient data in an
issue or email.

## License

Copyright (c) 2026 Wei Xu.

RepliAPA is released under the MIT License. See `LICENSE` for the complete
terms.
