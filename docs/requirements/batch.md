# External Batch Orchestration Requirements

**Requirement namespace:** `SRS-BAT-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-BAT-001:** The external batch runner MUST validate the selected manifest prefix, trace directory, reference, configuration, selected trace workload, identities, destinations, and cleanup targets before deleting any artifact.
- **SRS-BAT-002:** Batch preflight MUST reject invalid or duplicate selected IDs, missing selected traces, traces matching multiple selected samples, unsafe target types, path escapes, and symlinked traces or cleanup targets.
- **SRS-BAT-003:** Unless `--no-build` is selected, the release build MUST complete successfully before cleanup. In all modes the selected binary MUST be a regular file before cleanup.
- **SRS-BAT-004:** Cleanup MUST destructively remove only selected `results/<sample-id>/` directories and selected `$DNA_LOG_DIR/<sample-id>.log` files. It MUST preserve every unselected result directory and unrelated or legacy trace log.
- **SRS-BAT-005:** After cleanup, each selected trace MUST run through an isolated one-file core CLI invocation and each generated JSON MUST be placed without overwrite. Those helper invocations MUST write operational logs only inside their temporary workspaces. When every trace for a selected sample succeeds, the wrapper MUST run the sample command over the complete trace set, persist only `$DNA_LOG_DIR/<sample-id>.log`, and atomically publish the aggregate as `results/<sample-id>/<sample-id>.json`. The sample log MUST contain the nested trace-stage records produced by the authoritative read pipeline. The wrapper MUST NOT weaken the core CLI's no-overwrite semantics.
- **SRS-BAT-006:** If any trace for a sample fails, the batch runner MUST NOT publish that sample's aggregate JSON. Batch cleanup is not an all-workload transaction: a failure after cleanup MAY leave successful new results/logs from earlier traces or samples and MUST NOT claim to restore the removed prior selected artifacts; the final status MUST report trace and sample failures.
