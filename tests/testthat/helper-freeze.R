make_test_freeze <- function() {
  repliapa_create_freeze(
    config = list(model = "TEST", pseudocount = 0.5),
    pas_ids = c("p1", "p2"),
    primary_metrics = c("null_calibration", "replication"),
    role_registry = data.frame(
      registry_id = c("DEV", "REP2", "K562", "MPRA", "PATIENT", "NTC"),
      dataset = c("development", "validation", "external", "orthogonal",
                  "patient", "negative_control"),
      role = c("DEVELOPMENT", "FROZEN_INTERNAL_VALIDATION",
               "EXTERNAL_CELL_CONTEXT_VALIDATION",
               "OPTIONAL_ORTHOGONAL_VALIDATION",
               "METADATA_ONLY_APPLICATION_CANDIDATE",
               "NEGATIVE_CONTROL;NULL_CALIBRATION"),
      stringsAsFactors = FALSE
    ),
    frozen_on = "2026-08-11"
  )
}
