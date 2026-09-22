# Fuzzing

DNA treats ABIF/AB1 as untrusted binary input. The fuzz target exercises the
bounds-checked directory parser directly, without filesystem I/O.

Run locally on a supported Unix-like host:

```bash
cargo install cargo-fuzz --version 0.13.2 --locked
cargo +nightly fuzz run abif_parse -- -max_total_time=60
```

CI runs a short smoke campaign on pull requests that change Rust/parser/fuzz
code and a longer scheduled campaign. Any minimized input that causes a crash,
panic, timeout, or invariant violation must be retained as a regression fixture
before the underlying defect is considered fixed.
