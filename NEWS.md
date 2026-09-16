# RepliAPA 0.99.0 (v1.0.0-rc4, local candidate)

- Replaces the brief development README with a manuscript-grade GitHub user
  manual covering scope, installation, exact input/output contracts, workflow,
  parameters, diagnostics, reproducibility, limitations, reporting, data
  availability, citation and support without changing scientific inference.
- Finalizes package author, creator, maintainer, and copyright metadata for
  Wei Xu (徐炜; `xw222300@163.com`) while retaining the MIT license.
- This is a metadata-only release-candidate update: the rc4 API contract and
  all frozen scientific results are unchanged, and no public release occurred.
- Adds `repliapa_read_pas_counts()` for aligned matrix, Seurat,
  SingleCellExperiment, and SummarizedExperiment inputs.
- Requires an explicit PAS assay for containers and fails closed rather than
  silently reading a default RNA assay.
- Adds a complete executable vignette from input standardization and
  count-conserving pseudobulk through adequacy audit, fitting, calibration,
  frozen-replicate validation, all-PAS comparison, and verified result bundles.
- Adds machine-readable software citation metadata. Public repository and DOI
  fields remain absent until an owner-authorized external release.
- Extends the frozen API contract additively to 40 exports and 17 result/input
  classes; contract SHA256 is
  `29467b553b7df6ed26a21a8dac19eca3f223f6afac5f249d65eeeb3ed68addff`.

## v1.0.0-rc3 visualization changes

- Adds a publication-grade visual system without changing any estimand,
  fitted result, threshold, or frozen validation outcome.
- Exports `repliapa_theme()` and `repliapa_palette()` plus seven result-native
  ggplot2 functions for PAS effects, guide evidence, null calibration,
  biological-replicate validation, all-PAS/proximal-distal comparison,
  cross-context support, and APA-versus-gene-expression separation.
- Plotting is an optional layer: core inference retains its lightweight
  dependency set and every plot returns an editable ggplot object.
- Adds a worked visualization vignette covering biological interpretation,
  guide/replicate evidence, null calibration, representation comparisons,
  APA-versus-GE separation, and vector/raster publication export.
- Extends the fail-closed API contract additively to 39 exported functions;
  contract SHA256 is
  `a3fc9f773625446cfc1fa7e9d321a267b8c25b93446311296c2e1aada9f983f6`.

## v1.0.0-rc2 terminology changes

- Replaces the newly emitted abstention reason
  `TOO_MANY_EXACT_PERMUTATIONS` with
  `TOO_MANY_ADMISSIBLE_PERMUTATIONS`. The legacy value remains accepted by the
  contract for historical bundles but is no longer emitted.
- Describes inference as guide-aware finite-sample permutation calibration.
  Exhaustive enumeration is a computational property, not a claim of
  design-based exact randomization inference.
- This additive enum/terminology migration does not change estimands,
  statistics, p-values, thresholds, or frozen Phase 0 results.

- Promotes the reviewed v1 API/schema to a local release candidate using the
  R-valid pre-release version `0.99.0`; `v1.0.0-rc2` remains the explicit
  candidate label because hyphenated versions are invalid in `DESCRIPTION`.
- Contains no model, estimand, threshold, data-role, or scientific-result
  change relative to 0.11.1.
- External publication remained prohibited at this historical checkpoint until
  maintainer identity/contact metadata and a real version-control release tag
  were available. Maintainer metadata has since been finalized.

# RepliAPA 0.11.1

- Completes the v1 release-candidate API review without changing exported
  function signatures or valid-input scientific results.
- Rejects overlapping development/validation biological replicates and requires
  identical guide units in all-PAS versus proximal/distal comparisons.
- Fails early on global batch configuration errors instead of recording every
  target as biologically unsupported.
- Validates finite thresholds, exact integer resampling controls, target labels,
  replicate availability, gene-expression identifiers, and named fit routing.
- Allows arbitrary non-empty replicate suffixes in declared `effect_<replicate>`
  bundle extensions, including realistic underscore-containing labels.

# RepliAPA 0.11.0

- Freezes a fingerprinted v1-candidate contract for all exported function
  signatures, public S3 result components, bundle tables, enumerations, and
  compatibility rules.
- Adds `repliapa_api_contract()` for machine-readable inspection and
  `repliapa_verify_contract()` for fail-closed API and optional result-bundle
  verification without rerunning inference.
- Canonicalizes the zero-row target-effect schema so supported and unsupported
  result bundles expose the same required columns.
- Adds a full v0.1-v0.10 historical integration regression gate; no estimand,
  frozen threshold, or scientific result is changed.

# RepliAPA 0.10.0

- Adds validated externally instrumented compute tables to auditable result
  bundles; absent measurements still remain explicitly uninstrumented.
- Benchmarks four all-PAS/two-bin RepliAPA and DRIMSeq configurations in three
  independent R processes on one checksum-locked frozen 250-gene contrast.
- Separates model elapsed scope from whole-process kernel VmHWM peak RSS scope.
- Confirms an all-PAS compute advantage for RepliAPA while retaining the
  proximal/distal runtime trade-off where DRIMSeq is faster.

# RepliAPA 0.9.0

- Adds canonical target, PAS, guide, heterogeneity, unsupported-target, and
  tested-gene accounting tables without rerunning inference.
- Adds deterministic TSV result bundles with per-table dimensions, SHA256
  manifests, verification, and tamper detection.
- Result writing fails closed on non-empty directories and never removes
  unrelated user files, even when overwrite is requested.
- Adds optional verified-freeze and complete input-audit provenance. Runtime and
  peak memory remain explicitly `NOT_INSTRUMENTED` when no measurement exists.

# RepliAPA 0.8.0

- Adds a pre-fit pseudobulk data-contract audit without outcome access.
- Adds sample-unit qualification, guide/replicate support, and target eligibility
  with independent guides required in every biological replicate.
- Adds predeclared-development-only PAS-A0/A1/A2/A3 support tiers and explicit
  multi-PAS gene eligibility, including exact coordinate-conflict accounting.
- Reproduces the complete frozen HEK adequacy audit across 284 sample units,
  35,882 PAS, 12,617 genes, and all targets present in the processed object.

# RepliAPA 0.7.0

- Adds deterministic, count-conserving DRIMSeq input export from the identical
  guide-by-replicate pseudobulk contrast used by RepliAPA.
- Exports no hidden filter and records sample, guide, PAS, gene, and count
  accounting before any baseline-native filtering.
- Adds fair baseline comparison with separate method-native and matched tested-
  gene summaries; unsupported or untested genes remain explicit.
- Reproduces the frozen 250-gene NUDT21 DRIMSeq and RepliAPA results to machine
  precision through the new interfaces.

# RepliAPA 0.6.0

- Adds SHA256-fingerprinted analysis freeze objects containing model settings,
  fixed PAS IDs, primary metrics, data roles, and governance-source hashes.
- Adds fail-closed data-role/action authorization for development, frozen
  internal validation, external context, NTC, orthogonal, and metadata-only data.
- Adds a freeze-authorized external-context entry point and replaces the weak
  MPRA boolean gate with a verified freeze plus orthogonal registry identity.
- Adds payload mutation and governance-source tamper detection.

# RepliAPA 0.5.0

- Adds proximal-direction extraction from frozen all-PAS effect vectors.
- Adds post-freeze orthogonal proximal/distal count validation with within-
  replicate control pairing.
- Adds construct mapping, replicate adequacy, and unsupported accounting;
  orthogonal outcomes cannot be accessed through the API before model freeze.

# RepliAPA 0.4.0

- Adds independent cross-cell-context effect alignment on common PAS.
- Adds frozen `SHARED_EFFECT`, `OPPOSITE_EFFECT`,
  `CELL_LINE_SPECIFIC_EFFECT`, and `SUPPORTED_NEITHER_SIGNIFICANT` categories.
- Adds explicit target- and gene-level unsupported accounting so inadequate
  external guide support is not counted as biological non-replication.

# RepliAPA 0.3.0

- Adds a direction-aware proximal/distal collapse and matched two-bin fit.
- Adds explicit all-PAS versus two-bin tested-set comparison.
- Adds deterministic same-target guide reproducibility against a
  covariate-matched different-target null, with guide-pair and target bootstrap.

# RepliAPA 0.2.0

- Adds multi-target fitting with explicit unsupported-target accounting.
- Adds generated or predeclared guide-unit NTC null calibration.
- Adds equal-guide, per-cell RNA effects and APA/GE four-quadrant classification.

# RepliAPA 0.1.0

- Initial transparent implementation derived from the completed Phase 0 audit.
- Adds pseudobulk construction, exhaustive admissible guide permutation, PAS effect vectors,
  guide heterogeneity, and frozen-replicate reproducibility summaries.
