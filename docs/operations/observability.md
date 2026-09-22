# Observability

DNA uses typed failures and separate append-only operational logs to make command execution inspectable without changing deterministic scientific result contracts.

## Current signals

- typed stage-specific errors at the application boundary;
- append-only per-operation/sample logs;
- elapsed stage timing where recorded by the pipeline;
- deterministic result artifacts and source/config/reference identities used for reproducibility.

Scientific result JSON is not an operational log. Logging must not silently add dense identifying sequence/signal payloads.

When operational metrics or a continuously running service boundary is introduced, define explicit SLI/SLO documents rather than retrofitting service semantics onto the current CLI.
