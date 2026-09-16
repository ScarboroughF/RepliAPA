# Development status

## 2026-08-19 — manuscript-grade GitHub documentation

- Replaced the development README with a full English user manual aligned to
  software-paper reproducibility expectations.
- The manual distinguishes upstream PAS calling from RepliAPA inference,
  specifies statistical units and aligned input schemas, documents major
  parameters and result components, provides executable workflows, and states
  validation boundaries and reporting requirements.
- Stable repository and archive fields remain explicit `TO_BE_ADDED` items;
  no URL, release or DOI was invented or published.
- This documentation-only change does not modify the rc4 API contract or any
  frozen scientific result.

## 2026-08-19 — owner metadata finalized

- Package author, creator, maintainer, and copyright holder are Wei Xu
  (徐炜; `xw222300@163.com`).
- The existing OSI-approved MIT license is retained.
- No public repository, release tag, or archive DOI was created; external
  publication remains a separate owner-authorized action.
- This metadata-only change does not modify the rc4 API contract, estimands,
  statistics, thresholds, tests, or frozen scientific results.

## 2026-08-19 — v1.0.0-rc4 input and documentation contract

- Added one standardized, fail-closed reader for matrix, Seurat,
  SingleCellExperiment, and SummarizedExperiment PAS assays.
- Container input requires an explicit PAS assay name, preventing a default RNA
  assay from being misinterpreted as poly(A)-site counts.
- Added an executable end-to-end vignette covering input, cell-to-guide
  pseudobulk, audit, fitting, null calibration, frozen replicate validation,
  all-PAS/two-bin comparison, visualization, and verified bundles.
- Added package citation metadata without inventing a public DOI, repository,
  or author affiliation; owner identity was finalized subsequently.
- Extended the API contract to 40 exports and 17 classes under SHA256
  `29467b553b7df6ed26a21a8dac19eca3f223f6afac5f249d65eeeb3ed68addff`.
- This additive software-interface change does not modify estimands, statistics,
  thresholds, or frozen scientific results.

## 2026-08-19 — v1.0.0-rc3 visualization contract

- Added a coherent publication visual grammar and nine additive public
  visualization exports (theme, palette, and seven native result plots).
- Visualizations consume existing frozen result objects and do not rerun,
  refit, tune, filter, or reinterpret scientific inference.
- Extended the candidate API contract from 30 to 39 exports under SHA256
  `a3fc9f773625446cfc1fa7e9d321a267b8c25b93446311296c2e1aada9f983f6`.
- Added visualization unit tests covering native fit, guide, null,
  reproducibility, all-PAS, context, and APA/GE outputs.

## 2026-08-13 — v1.0.0-rc2 terminology contract

- Replaced the emitted abstention machine code
  `TOO_MANY_EXACT_PERMUTATIONS` with
  `TOO_MANY_ADMISSIBLE_PERMUTATIONS` to prevent exhaustive computation from
  being misread as design-based exact randomization inference.
- The rc2 contract accepts the historical code when validating older bundles,
  while new fits emit only the new code. This is an explicitly disclosed,
  additive enum migration and does not change scientific inference.
- Candidate contract SHA256 is
  `cc8286e79df46e5887a769c8024ce93219607083cb03d83651ec7807c31f0906`.
- Package regression tests and formal check status must be rerun before rc2 can
  be considered release-ready.

## 2026-08-13 — local v1.0.0-rc1 candidate

- Promoted the reviewed API/schema to R-valid package version `0.99.0`, with
  explicit candidate label `v1.0.0-rc1` and `LOCAL_ONLY` release scope.
- R does not accept `1.0.0-rc1` in `DESCRIPTION`; no Git executable is present,
  so an external checksum-protected marker is used and is not called a Git tag.
- Contract, full bundle, and unsupported bundle checks pass 45/45, 77/77, and
  61/61. Historical v0.1-v0.10 regression remains 10/10 PASS and scientific
  inference remains unchanged.
- `R CMD check --no-manual RepliAPA_0.99.0.tar.gz`: **Status: OK**, with zero
  errors, warnings, or notes.
- The external release manifest records the final tarball SHA256 to avoid a
  self-referential package artifact. Candidate contract SHA256 is
  `849dc36688762decf55acbb9b86acf5ff60d4d089f64562a5e435f103786a02f`.
- No external publication occurred. At that historical checkpoint, placeholder
  maintainer metadata and absent Git tagging remained external-release blockers.

## 2026-08-13 — v0.11.1 release-candidate API review

- Reviewed all 30 exported functions, 16 public result classes, 15 canonical
  bundle tables, documentation topics, defaults, failure behavior, and named
  integration routing without changing a public function signature.
- Added mandatory disjoint-replicate enforcement to biological-replicate
  reproducibility and identical-guide enforcement to the all-PAS/two-bin
  comparison, preventing invalid validation claims.
- Global batch configuration errors now fail before target iteration rather
  than being mislabeled as biological target non-support.
- Added explicit finite-number, exact-integer, target, replicate, threshold,
  gene-expression identity, and named-fit routing guards. Fractional resampling
  counts are no longer silently truncated.
- Corrected the contract extension rule so valid `effect_<replicate>` columns
  accept realistic labels containing underscores or other non-empty suffixes.
- Re-ran all v0.1-v0.10 real-data integrations at full scale: 10/10 PASS. The
  v1 candidate full/unsupported bundle checks are 77/77 and 61/61 PASS;
  `scientific_inference_changed = FALSE`.
- Final candidate contract SHA256 is
  `76a7347447b17dce807b61f07fa8f225abfdce6190efb9def808b56a2e87f561`.
- `R CMD check --no-manual RepliAPA_0.11.1.tar.gz`: **Status: OK** with zero
  errors, warnings or notes. The external release manifest records the final
  source-tarball checksum to avoid a self-referential package artifact.

## 2026-08-12 — v0.11 v1-candidate API/schema freeze

- Added a machine-readable contract for 30 exported function signatures, 16
  public result classes, 15 canonical bundle tables, declared extension
  patterns, storage types, enums, and compatibility policy.
- The contract is fail-closed under SHA256 fingerprint
  `17ab1302...0c0a60`; payload edits, signature changes, missing/undeclared
  columns, type changes, and undeclared PAS-effect extensions are detected.
- Added contract inspection and verification APIs. Verification never reruns
  inference and can validate both full 15-table bundles and zero-row
  unsupported-target bundles.
- Re-ran every historical integration from v0.1 through v0.10 without reducing
  tested genes or resampling. All 10 milestones passed, including 77,000 guide-
  gene pairs, 3,735 cross-context pairs, 746 MPRA constructs, 35,882 PAS, the
  DRIMSeq comparison, governance guards, result bundles, and compute provenance.
- All scientific differences remain at their prior machine-precision bounds;
  `scientific_inference_changed = FALSE`.
- `R CMD check --no-manual RepliAPA_0.11.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes. Source SHA256 is
  `5ec1849d9e6bcd3d29eeee2865a1177ce64cfe6c8d8deb3f20c7ecdb5788df32`.

## 2026-08-12 — v0.10 compute-cost milestone

- Added validated measured-compute provenance to canonical result bundles.
  Runtime and memory statuses change only when a complete external benchmark
  table supplies positive measurements, input SHA256, and explicit scopes.
- Prepared one checksum-locked NUDT21 benchmark input containing the identical
  1,082 PAS, 250 genes, 26 guide-by-replicate units, 3 target guides, and 10 NTC
  guides used by the frozen method comparison.
- Ran four methods in three sequential independent R processes each. Model
  elapsed time excludes input loading; peak RSS is the Linux kernel VmHWM for
  the whole isolated process and includes package/input residency.
- Median all-PAS elapsed/peak RSS were 8.583 s/195.1 MB for RepliAPA and
  13.804 s/606.5 MB for DRIMSeq. RepliAPA therefore did not buy its inferential
  advantage with higher compute cost on this frozen contrast.
- For proximal/distal, DRIMSeq was faster (8.235 s vs 11.207 s), while RepliAPA
  used lower process peak RSS (192.0 MB vs 594.5 MB). This runtime trade-off is
  retained rather than hidden.
- All three repeats per method have one deterministic result digest. RepliAPA
  all-PAS and both Phase 0 DRIMSeq methods reproduce frozen P/q/tested status to
  machine precision; package proximal/distal is marked determinism-only because
  no Phase 0 golden test table exists for that exact engine/representation pair.
- `R CMD check --no-manual RepliAPA_0.10.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.9 auditable-result-bundle milestone

- Added canonical target, PAS, long guide-level PAS, heterogeneity, target-
  support, and tested-gene accounting tables for single or multi-target fits.
- Unsupported targets, filtered genes, abstentions, and native tested fractions
  remain explicit. Bundle construction never reruns or changes inference.
- Added deterministic TSV writing, table dimensions, SHA256 manifests, bundle
  verification, tamper detection, and fail-closed overwrite behavior that
  refuses directories containing unrelated user files.
- The formal frozen NUDT21 bundle contains 14 tables: 250 target-effect rows,
  1,082 PAS-effect rows, 6,492 guide-level PAS-effect rows, and 500
  heterogeneity rows. All 250 genes are tested with zero filtering/abstention.
- The result bundle carries the verified Phase 0 freeze fingerprint and complete
  HEK input audit. Statistics/P/q/effect norms reproduce the reference within
  `1.1e-15`; coordinate shifts differ by at most `1.1e-11`.
- Runtime and peak memory are truthfully recorded as `NOT_INSTRUMENTED`; no
  retrospective values are fabricated.
- `R CMD check --no-manual RepliAPA_0.9.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.8 input-contract milestone

- Added an outcome-free audit for aligned PAS pseudobulk counts, sample-unit
  metadata, PAS-to-gene mapping, and the no-cell-as-replicate contract.
- Added transparent sample-unit qualification and guide/target support. Primary
  eligibility requires at least two qualified independent guides in every
  biological replicate; unsupported targets remain explicit.
- Added PAS-A0/A1/A2/A3 tiers calculated only from predeclared development
  replicates and qualified units, plus multi-PAS gene eligibility and exact
  coordinate-conflict accounting.
- The full frozen HEK regression exactly reproduces all support tiers and
  eligibility calls: 270/284 qualified sample units, 33/36 eligible non-NTC
  targets present in the processed object, 35,401/35,882 A2/A3 PAS, and
  8,473/12,617 primary multi-PAS genes. PAS, gene, and target mismatches are 0.
- One target in the later registry but absent from the processed object remains
  reference-only and is not fabricated by the input audit. Outcome access is
  recorded as `FALSE`.
- `R CMD check --no-manual RepliAPA_0.8.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.7 baseline-interoperability milestone

- Added deterministic DRIMSeq-standard count and sample export from the exact
  target/NTC guide-by-replicate contrast used by RepliAPA.
- The export performs no baseline filtering, preserves a predeclared PAS order,
  retains complete guide identities across replicates, and reports count,
  sample, PAS, gene, and exclusion accounting.
- Added baseline result alignment with separate method-native and common tested-
  gene summaries, discovery accounting, rank concordance, and overlap.
- The frozen NUDT21 regression exports 1,082 PAS from 250 genes across 26 sample
  units (3 target and 10 NTC guides in two replicates) with exact count
  conservation. DRIMSeq accepts the exported object directly.
- All 250 DRIMSeq and RepliAPA gene results reproduce the Phase 0 references;
  maximum P/q difference is below `5.6e-16`. One DRIMSeq-untested gene remains
  explicit, giving 249 genes in the matched comparison.
- `R CMD check --no-manual RepliAPA_0.7.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.6 executable-governance milestone

- Added an analysis freeze object containing model configuration, 1,082 frozen
  calibration PAS IDs, primary metrics, eight Phase 0 data-role records, source
  SHA256 values, and a deterministic payload fingerprint.
- Added freeze verification, payload mutation detection, optional live source-
  file verification, and fail-closed role/action authorization.
- HEK replicate 2 and K562 allow frozen application but prohibit model/filter/
  metric changes; MPRA requires its orthogonal role; patient candidates remain
  metadata-only; NTC permits guide-unit null calibration but prohibits splitting.
- Added a freeze-authorized external-context entry point. Orthogonal validation
  now requires the verified freeze and registry identity rather than a boolean.
- Phase 0 governance regression reproduces the recorded model-config SHA256
  `73fedc22...ffeafc` and primary-metrics SHA256 `78529f43...20ac1`; all 11
  role, source-integrity, and mutation guards pass.
- `R CMD check --no-manual RepliAPA_0.6.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.5 orthogonal-validation milestone

- Added proximal-direction extraction from frozen all-PAS effects using the
  direction-aware two-bin membership.
- Added an explicitly post-freeze orthogonal validation API; in v0.6 its
  original boolean gate was strengthened to a verified freeze contract.
- Added within-replicate control pairing, row-to-construct aggregation,
  replicate support, mapping, untested, and absent-target accounting.
- The full MPRA regression processes 366,780 raw rows and reproduces all 38
  frozen mapped constructs across 18 genes with zero direction mismatches;
  maximum construct-effect difference is below `4.8e-16` and maximum summary
  difference is below `3.4e-16`.
- One previously coarsely unmapped construct lacking a second MPRA replicate is
  now explicitly orthogonal unsupported; mapped NUDT21/CSTF3 conclusions are
  unchanged.
- `R CMD check --no-manual RepliAPA_0.5.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.4 external-context milestone

- Added common-PAS alignment of independently fitted primary and external
  contexts without external-outcome tuning.
- Added frozen shared, opposite, cell-line-specific, and neither-significant
  gene-level categories.
- Added explicit `EXTERNAL_UNSUPPORTED`, `PRIMARY_UNSUPPORTED`, and
  `BOTH_UNSUPPORTED` target accounting plus common-site gene abstention.
- HEK/K562 regression reproduces all 3,735 frozen target-by-gene categories for
  15 adequately supported targets with zero category mismatches; maximum
  effect-vector cosine difference is below `6.7e-16`.
- All 18 HEK-eligible targets lacking K562 support remain explicitly external
  unsupported and are not counted as biological failures.
- `R CMD check --no-manual RepliAPA_0.4.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.3 added-value and guide-validation milestone

- Added direction-aware PAS-to-proximal/distal collapse with exact count
  conservation and explicit single-PAS/invalid-strand accounting.
- Added a two-bin fit using the identical guide permutation engine as all-PAS,
  plus matched tested-set comparison between the representations.
- Added lossless collection of PAS-level guide effects and deterministic
  same-target guide comparison against covariate-matched different-target pairs.
- Matching stays within biological replicate, requires distinct null targets,
  and bootstrap resamples guide pairs or targets rather than cells.
- HEK NUDT21 regression reproduces all 250 frozen distal-bin effect summaries;
  maximum absolute difference is below `4.8e-16`.
- The complete frozen guide validation is reproduced: 308 guide pairs and
  77,000 pair-by-gene rows, with similarity differences below `6.7e-16`,
  matching-distance differences below `4.9e-15`, and bootstrap-CI differences
  below `4.5e-16`.

## 2026-08-12 — v0.2 audit-workflow milestone

- Added multi-target fitting with explicit fitted/unsupported accounting.
- Added generated and predeclared whole-guide NTC pseudo-target calibration;
  guide identity and biological-replicate blocks remain intact.
- Added equal-guide, cell-count-normalized RNA effects and explicit
  `APA_ONLY`, `GE_ONLY`, `APA_AND_GE`, and `NEITHER` classification.
- Unit tests cover batch accounting, declared NTC contrasts, deterministic
  contrast generation, caller RNG preservation, and APA/GE separation.
- HEK regression matches the frozen Phase 0 outputs for three NTC contrasts
  (750 gene-tests): maximum statistic/P/q differences are below `6e-16`.
- HEK NUDT21 gene-expression regression matches all 250 frozen effects with a
  maximum absolute log2FC difference below `2.3e-15`.
- The original NUDT21 all-PAS migration regression remains unchanged.
- `R CMD check --no-manual RepliAPA_0.2.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.

## 2026-08-12 — v0.1 core milestone

- Human continuation authorization received after Phase 0 RAPA-E3.
- Six public APIs implemented and documented.
- Unit tests cover count and sample-unit contracts, count conservation, exact
  guide permutations, compositional effects, evidence retention, filtering,
  abstention, determinism and replicate reproducibility.
- `R CMD check --no-manual RepliAPA_0.1.0.tar.gz`: **Status: OK** with zero
  errors, warnings or notes.
- HEK NUDT21 250-gene migration regression matches the frozen Phase 0 outputs:
  maximum statistic/p/q/effect-norm differences are below `6e-16`; maximum
  coordinate-shift difference is below `6e-12`.
This milestone does not authorize or implement TMB, deep models, patient raw
data analysis, a pkgdown site, or Bioconductor submission.
