# Code Quality

DNA prefers one clear implementation with explicit ownership and no transitional compatibility layer unless the current specification requires one.

Required engineering properties include:

- formatted, warning-free Rust and Python tooling;
- typed failures rather than hidden fallback behavior;
- no first-party unsafe Rust without a new explicit decision and targeted validation;
- no deprecated/legacy compatibility scaffolding or diagnostic suppression used to preserve obsolete paths;
- deterministic behavior where required by the domain;
- bounded resource use for untrusted inputs;
- focused modules with clear dependency direction;
- documentation and tests synchronized with changed behavior.

## Rust lint policy

`Cargo.toml` `[lints]` is the single owner of compiler and Clippy policy; crate
roots carry no lint attributes. CI runs Clippy with `-D warnings`, so every
`warn` level is a gate.

- rustc: `unsafe_code` is forbidden, `deprecated` is denied, and
  `missing_docs`, `unreachable_pub`, `unused_qualifications`, and
  `rust_2018_idioms` are enforced.
- Clippy: `unwrap_used` and `expect_used` are denied and the `pedantic` group is
  enforced, except three lints allowed crate-wide with a recorded rationale:
  `cast_precision_loss` (counts are bounded by the input limits in
  `src/config/defaults.rs`), `float_cmp` (byte-deterministic outputs make exact
  float assertions intentional), and `too_many_lines` (stage functions mirror
  one documented stage record).
- A site that is provably safe but still trips a lint uses a targeted
  `#[expect(lint, reason = "...")]` that states the proof. `expect` fails the
  build if the exception becomes unnecessary. Suppressing dead, unused,
  unreachable, or deprecated code is never permitted (INV-RUST-004).

The executable quality gates are defined by [CI/CD](ci-cd.md) and repository policy validators.
