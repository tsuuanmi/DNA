# `src/report/atomic.rs`

## Purpose

Publishes one completed JSON file atomically without overwriting an existing
target.

## Responsibilities

- Create the parent directory if needed, reserve a collision-free temporary
  sibling through bounded PID/counter attempts, write and sync it, then hard-link
  it to the target and remove the temporary.
- Sync the parent directory after publication.
- Clean up the temporary file on failure.

## Non-responsibilities

No serialization, format logic, or algorithm execution.

## Key types and functions

- `publish(path, bytes) -> Result<()>`: the module entry point, re-exported from
  `mod.rs`.
- `create_temporary`: bounded collision-resistant sibling reservation.
- `TemporaryFile`: a RAII guard that removes the temporary file unless published.

## Invariants and errors

- The target must not already exist; otherwise `Error::Path`.
- The parent directory is created as needed; a failure to create or open it
  returns `Error::Output`.
- File operations return `Error::Output` on failure.
- The temporary file is removed on any failure path.

## Dependencies

- `error` for `Error`/`Result`.

## Biological semantics

None; this module is purely about safe output publication.

## Tests

- `skips_a_stale_temporary_name` verifies bounded collision handling.
- The integration test `refuses_to_overwrite_completed_output` verifies that an
  existing target is not overwritten.

## Status

Implemented.
