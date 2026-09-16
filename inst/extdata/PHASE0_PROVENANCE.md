# Phase 0 provenance

The statistical core was frozen before HEK replicate-2 target outcomes in
`RepliAPA_phase0/00_governance/MODEL_CONFIG_FREEZE.yaml` with SHA256
`73fedc22a677b7be2f067b874e799a1540abcd55d29b13490928c2bde2ffeafc`.

Phase 0 concluded at RAPA-E3 and `PACKAGE_DEVELOPMENT_GO`. This package does
not bundle the public CPA-Perturb-seq data or its processed objects.

Package migration regressions cover the core target fit (v0.1), NTC and RNA
effects (v0.2), and direction-aware PAS bins plus all 77,000 frozen matched
guide-pair/gene comparisons (v0.3).
The v0.4 regression additionally covers all 3,735 frozen HEK/K562 common-site
target-by-gene classifications and the full target-support universe.
The v0.5 regression covers all 366,780 MPRA raw rows, 746 construct-target
records, and the 38 mapped NUDT21/CSTF3 direction comparisons.
The v0.6 regression verifies the original model/metric hashes, creates an
executable eight-record role registry, and passes all 11 frozen-action and
tamper guards.
The v0.7 regression passes the frozen 1,082-PAS NUDT21 contrast through the
DRIMSeq-standard export, reproduces all 250 DRIMSeq and RepliAPA gene results
to machine precision, and retains the 249-gene matched tested set separately
from each method's native accounting.
The v0.8 regression reproduces the frozen HEK input contract and adequacy audit
for all 284 sample units, 35,882 PAS, 12,617 genes, and every target actually
present in the processed object, without reading target outcomes.
The v0.9 regression writes and verifies a 14-table NUDT21 result bundle linked
to the Phase 0 freeze and complete input audit. All fitted statistics remain
unchanged, and absent compute instrumentation is reported rather than inferred.
The v0.10 audit adds three isolated-process repetitions for four method/
representation pairs on a checksum-locked frozen input. It measures model-scope
elapsed time and kernel VmHWM peak RSS, while preserving the original Phase 0
whole-pipeline rows as historically uninstrumented.
