# Analysis Failure Investigation

Use this playbook when a DNA command fails and the cause is not already known.

## 1. Preserve evidence

- keep the terminal typed error and corresponding append-only log;
- preserve input/reference/config identities without copying sensitive biological payloads into issue text;
- do not rerun with relaxed validation or compatibility fallbacks.

## 2. Classify the failing boundary

Determine whether the failure is in:

- input/path/config validation;
- ABIF or FASTA decoding;
- basecalling/signal/QC;
- alignment or variant processing;
- sample aggregation;
- serialization/publication;
- operational logging/filesystem behavior.

## 3. Reproduce with the narrowest safe case

Use synthetic/non-sensitive data where possible. Confirm the exact source revision, configuration, reference identity, and command.

## 4. Compare against canonical intent

Read the relevant [requirements](../../requirements/README.md), [design](../../design/README.md), [reference](../../reference/README.md), and source README. A disagreement between implementation and normative docs is a defect to reconcile.

## 5. Escalate

If the issue is a scientific disagreement, follow [validation strategy](../../validation/strategy.md). If it is a security/data concern, follow [security](../../security/README.md) and [data governance](../../governance/data.md).
