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

The executable quality gates are defined by [CI/CD](ci-cd.md) and repository policy validators.
