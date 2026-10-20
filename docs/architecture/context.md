---
decided-by: [ADR-0065, ADR-0066, ADR-0075]
---

# System Context

DNA is a local scientific Rust library and CLI for DNA analysis.

The current production capability processes Sanger sequencing traces stored as
ABIF. The platform architecture is intentionally broader than that implementation
and is not restricted to mitochondrial DNA.

## Context diagram

```mermaid
flowchart LR
    subgraph inputs["inputs (untrusted)"]
        ab1["AB1 traces"]
        ref["reference FASTA"]
        cfg["strict TOML + target profile"]
        docs["JSON documents:<br/>read evidence, consensus, variants, run signal"]
        cf["reviewed consensus FASTA"]
    end

    dna["DNA<br/>library and CLI (dna)<br/>deterministic, local, no network"]

    subgraph outputs["outputs"]
        json["one versioned JSON result per command<br/>(+ consensus FASTA)"]
        logs["append-only operation logs"]
        api["typed Rust results<br/>(public capabilities)"]
    end

    ab1 --> dna
    ref --> dna
    cfg --> dna
    docs --> dna
    cf --> dna
    dna --> json
    dna --> logs
    dna --> api

    json --> review["laboratory review tool"]
    json --> pipelines["downstream pipelines<br/>orchestrate, convert, compare"]
    eval["dna-eval<br/>runs DNA over a panel,<br/>compares with reviewer truth"] --> dna
    truth["reviewer truth"] --> eval
```

## Platform scope

DNA separates four orthogonal concerns:

- **biological target** — for example mitochondrial DNA, nuclear/genomic DNA, or
  targeted loci/panels;
- **sequencing modality** — currently Sanger traces and reviewed consensus
  sequences, with NGS as a future modality;
- **external format** — for example ABIF, FASTA/FASTQ, BAM/CRAM, or VCF/BCF;
- **analysis capability** — for example base calling, alignment/mapping, variant
  analysis, SNP/genotyping analysis, haplogroup, or nomenclature.

SNP is therefore an analysis target/variant class rather than another sequencing
modality.

## Current production inputs

- Sanger sequencing traces encoded as ABIF (commonly named `.ab1`);
- reviewed consensus sequences in FASTA, one record per read (`call`);
- recorded contract documents that replay a later step without the traces:
  `dna.read_evidence/v1` (`analyze --evidence`), `dna.consensus/v1` (`call`), and
  `dna.signal/v1` (`review --signal`);
- `dna.variants/v1` documents for post-calling notation (`notation`);
- one short FASTA reference for reference-guided operations;
- strict TOML configuration and the target profile it names;
- explicit sample identifiers and trace paths for sample aggregation.

Future formats or modalities are not production behavior until their own
requirements, adapters, validation, and public contracts are implemented.

## External outputs

- one versioned command-specific JSON result per successful core invocation;
- typed Rust results for public library capabilities;
- separate append-only operational logs;
- downstream pipelines orchestrate runs and group, convert, and compare outputs
  (ADR-0065, ADR-0066); the `dna` binary and library do none of that. The
  evaluation harness `dna-eval`, a separate workspace binary, runs DNA over
  test panels and compares with reviewer truth as validation evidence
  (ADR-0075).

## Trust boundary

Input bytes, paths, filenames, manifests, configuration, future HTS records, and
external scientific data are untrusted. Sensitive biological data is governed by
[[governance/data|data policy]] and [security](../security/README.md).

The core is not a network service and currently has no deployment topology. A
deployment document should be added only if a real deployed runtime boundary is
introduced.
