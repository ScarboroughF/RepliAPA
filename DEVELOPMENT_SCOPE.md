# v0.11 development contract

Public API:

1. `repliapa_read_pas_counts()` standardizes matrix, Seurat,
   SingleCellExperiment, and SummarizedExperiment PAS inputs without guessing
   an assay; `repliapa_validate_counts()` validates nonnegative integer counts.
2. `repliapa_pseudobulk()` aggregates cells into cell-line × replicate × guide.
3. `repliapa_validate_pseudobulk()` enforces one sample unit per guide/replicate.
4. `repliapa_fit()` estimates all-PAS effects and finite-sample,
   guide-permutation-calibrated tests.
5. `repliapa_reproducibility()` compares separately fitted frozen replicates.
6. `repliapa_example()` supplies a deterministic, non-biological example.
7. `repliapa_fit_targets()` batches targets without hiding unsupported targets.
8. `repliapa_ntc_calibration()` runs generated or predeclared guide-unit nulls.
9. `repliapa_gene_expression_effect()` and `repliapa_classify_effects()` keep
   total expression separate from APA composition.
10. `repliapa_collapse_proximal_distal()`, `repliapa_fit_proximal_distal()`, and
    `repliapa_compare_allpas_pd()` provide a controlled all-PAS/two-bin comparison.
11. `repliapa_collect_guide_effects()` and
    `repliapa_guide_reproducibility()` retain long PAS effects and compare
    same-target guides with covariate-matched different-target guides.
12. `repliapa_compare_contexts()` aligns independently fitted contexts on
    common PAS, classifies shared/opposite/context-specific effects, and keeps
    target- and gene-level unsupported cases explicit.
13. `repliapa_proximal_directions()` derives frozen direction estimates and
    `repliapa_orthogonal_validation()` performs post-freeze construct-level
    direction checks with within-replicate controls and explicit mapping.
14. `repliapa_create_freeze()`, `repliapa_verify_freeze()`, and
    `repliapa_authorize()` fingerprint the analysis contract and enforce data
    roles; `repliapa_apply_external()` records freeze provenance on external use.
15. `repliapa_drimseq_input()` exports the identical unfiltered sample-unit
    contrast to DRIMSeq, while `repliapa_compare_baseline()` reports both
    method-native and matched tested-gene comparisons.
16. `repliapa_audit_input()` performs outcome-free sample-unit, guide/target,
    PAS-tier, and multi-PAS gene eligibility audits before fitting.
17. `repliapa_result_bundle()`, `repliapa_write_bundle()`, and
    `repliapa_verify_bundle()` standardize result tables, preserve tested/
    filtered/abstained accounting, and protect deliverables with SHA256.
18. `repliapa_result_bundle(..., compute_benchmark = )` validates measured
    elapsed/peak-RSS provenance and carries it into the checksum-protected
    bundle without estimating missing compute values.
19. `repliapa_api_contract()` returns the fingerprinted v1-candidate API and
    result-schema contract; `repliapa_verify_contract()` fails closed on
    altered signatures, S3 registration, table names, columns, storage types,
    or undeclared extensions.

Compatibility boundary: rc4 freezes the current 40 exported signatures, 17
public result-class component sets, 15 canonical result-bundle tables, and
published enum meanings. Maintenance in this line may fix implementation bugs
but may not silently change the scientific estimand, Phase 0 thresholds, data
roles, or historical result values.

Deferred: PAS detection, MAAPER input construction, internal DRIMSeq fitting,
TMB mixed models, patient applications, pathway analysis, neural models,
websites and Bioconductor submission.
