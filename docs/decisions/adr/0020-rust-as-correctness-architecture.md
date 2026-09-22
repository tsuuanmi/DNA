# ADR-0020: Use Rust as correctness architecture, not only as an implementation language

## Status

Accepted

## Context

DNA previously had useful reference points at both ends of the implementation spectrum:

- exploratory Sanger/mtDNA work implemented in Python, where iteration is fast and the scientific ecosystem is broad;
- Tracy, a mature C++ Sanger chromatogram tool with strong performance, broad functionality, established distribution, and years of practical use.

Both are valuable. DNA does not use Rust because Python or C++ are incapable of production software.

The reason is architectural.

> As systems become too complex for humans to reliably maintain implicit invariants, software engineering shifts those invariants into forms machines can verify. Rust is growing not merely because it is fast or memory-safe, but because its design moves correctness from developer discipline into the programming model itself.

A Sanger/mtDNA pipeline accumulates many interacting invariants:

```text
ABIF offsets must be in bounds
channel order must remain canonical
PLOC coordinates must remain valid
call indexes and sample indexes are different units
reference coordinates and trace coordinates are different units
reverse-strand mappings must preserve original call identity
linear and circular references have different normalization rules
raw evidence must not be silently mutated by downstream processing
an unresolved observation is different from an absent observation
read-level evidence is different from sample-level evidence
a mixed signal is not automatically heteroplasmy
published output must correspond to one complete validated analysis
```

If these rules exist mainly in comments, conventions, reviewer memory, or tests, the cost of safe change rises with system complexity.

DNA therefore treats Rust as a way to move as many mechanically checkable invariants as practical into types, ownership, constructors, exhaustive matching, module boundaries, and compiler-enforced APIs.

## Decision

DNA will design its Rust implementation around **invariant migration**:

> Every invariant that can reasonably be represented in the programming model should be moved out of human memory and into a machine-checkable form.

This is not limited to memory safety.

### 1. Invalid states should be difficult or impossible to represent

Prefer validated domain types and constructors over primitive values with comments.

Examples of the intended direction include:

```text
CallIndex         rather than an unqualified usize
SampleIndex       rather than an unqualified usize
ReferencePosition rather than an unqualified usize
LinearReference | CircularReference
Forward | Reverse
CanonicalBase | UnresolvedBase
ValidatedTrace
ValidatedReference
NormalizedVariant
```

Not every primitive needs a wrapper. Newtypes and enums are justified when they prevent coordinate, strand, topology, state, or semantic confusion.

A type should represent a scientific fact only after the required validation has occurred.

### 2. Pipeline stages should exchange validated typed states

A later stage should not be able to accidentally consume an object that has not passed the earlier stage's invariants.

Conceptually:

```text
AB1 bytes
   |
   v
DecodedTrace
   |
   v
BaseCalledRead
   |
   v
QualityControlledRead
   |
   +--> BasecallResult
   |
   v
AlignedRead
   |
   v
VariantEvidence
   |
   v
AnalysisResult
```

This does not require an elaborate generic typestate framework. Ordinary Rust structs, private fields, module visibility, and validated constructors are preferred when they express the boundary clearly.

The goal is not type-level cleverness. The goal is to make incorrect stage composition harder than correct composition.

### 3. Source evidence is immutable by default

Decoded chromatogram channels are source evidence.

Downstream processing should borrow immutable evidence and return a new typed value:

```rust
fn process(trace: &DecodedTrace) -> Result<ProcessedTrace, DNAError>
```

rather than using broad in-place mutation.

Baseline correction, smoothing, denoising, locus refinement, alignment, and consensus must not destroy the evidence from which they were derived.

Rust's ownership and borrowing model makes this architectural preference enforceable at API boundaries rather than merely conventional.

### 4. Failure is explicit data

Expected failures from external input or scientific eligibility rules should use typed `Result`/`Option`/enums rather than:

- sentinel integers;
- magic strings;
- null pointers;
- partially initialized objects;
- implicit fallback;
- printing an error and continuing.

Untrusted ABIF data must never rely on panic behavior as validation.

An error type should carry the stage and reason necessary for deterministic handling without exposing private scientific payloads.

### 5. Exhaustive state handling is preferred

Biological and pipeline states should use enums when the alternatives form a closed set.

Examples:

```text
Forward | Reverse
Linear | Circular
SNV | Insertion | Deletion
Uncovered | CoveredNoCall | SingleRead | Bidirectional | Discordant
Primary | Ambiguous | Unresolved
```

When a new state is added, exhaustive matching should force affected code to make an explicit decision.

This is preferable to stringly typed states that silently flow through old code.

### 6. Resource and coordinate safety are part of correctness

Memory safety alone is not enough.

All externally derived lengths, counts, offsets, products, matrix dimensions, and allocations must remain bounded and checked before use.

The design should prefer APIs that make unit conversions explicit and preserve mappings rather than repeatedly reconstructing coordinates from raw integers.

### 7. Deterministic scientific transformations are first-class

Scientific stages should be pure or observational whenever practical:

```text
validated input + explicit configuration
            |
            v
      typed deterministic output
```

Filesystem publication, timestamps, logging, and operating-system state belong at orchestration boundaries.

This separation enables exact regression testing, property testing, fuzzing, mutation testing, and reproducible reasoning about scientific changes.

### 8. Safe Rust is the default trust boundary

First-party production code forbids `unsafe`.

If a future requirement genuinely needs unsafe Rust or FFI, that boundary requires an explicit ADR describing:

- why safe Rust is insufficient;
- the exact safety invariant the compiler cannot check;
- the smallest possible unsafe surface;
- tests and Miri/sanitizer strategy where applicable;
- dependency and platform implications.

The presence of an unsafe dependency is not automatically a defect, but it increases the review and supply-chain burden.

## Why this is preferable to the previous Python direction

Python remains highly valuable for:

- exploratory scientific analysis;
- notebooks and visualization;
- rapid experiments;
- corpus analysis;
- one-off data transformation;
- model training and statistics.

DNA already uses Python appropriately for external development/validation orchestration.

The distinction is that Python's runtime does not enforce function or variable type annotations. Static type checkers can provide strong assistance, and runtime validation libraries can add contracts, but those protections are layered on rather than fundamental to execution semantics.

For a production scientific core, this tends to leave more invariants distributed across:

```text
type checker configuration
runtime validation
tests
conventions
dataframe/array shape assumptions
environment/package management
developer discipline
```

Python can absolutely be made production-grade. For DNA's small deterministic scientific core, however, Rust lets more of those invariants become mandatory properties of the compiled program while still providing native performance.

Python remains a good companion language around the Rust core rather than the authority for binary parsing, coordinate-sensitive scientific state, or release output.

## Why this is preferable to simply adopting the Tracy C++ architecture

Tracy is an important scientific and engineering reference. It has mature Sanger functionality and substantially more deployment history than DNA.

DNA should learn from Tracy's algorithms and product capabilities without copying all of its implementation constraints.

C++ can express strong types, RAII, immutability, bounds-aware containers, expected-like error values, sanitizers, static analysis, and excellent tests. A disciplined modern C++ project can reach very high assurance.

The difference is the default enforcement model.

C++ permits pointer arithmetic, unchecked indexing, invalid references, use-after-free, data races, implicit conversions, and other forms of undefined behavior unless project discipline and additional tooling prevent them. Tracy's ABIF implementation, for example, necessarily performs direct integer offset and buffer manipulation and reports several input failures through conventional return values and stderr.

Rust moves a larger subset of these obligations into the compiler and standard programming model.

For DNA specifically, Rust also gives one coherent toolchain for:

```text
build
dependency resolution
lockfile
formatting
linting
testing
documentation
cross compilation
fuzz integration
release metadata
```

The current DNA dependency graph is also Rust-native and does not require Tracy's Boost/HTSlib/SDSL-style system dependency stack. That simplifies building and distributing a small standalone command-line scientific tool.

This is not a claim that Rust is universally superior to C++. It is a claim that DNA values **enforced invariants and a small auditable deployment surface** more than access to the broadest existing C++ bioinformatics ecosystem.

## What Rust does not prove

Choosing Rust does not establish scientific correctness.

The compiler cannot determine whether:

- a peak threshold is biologically justified;
- PLOC is the best locus authority;
- an alignment scoring model is appropriate;
- a mixed peak is artifact, contamination, mixture, or heteroplasmy;
- a quality metric is calibrated;
- a reported mtDNA difference is clinically meaningful;
- a validation corpus is representative.

Those are empirical and methodological questions.

Rust also does not eliminate:

- logic bugs;
- incorrect formulas;
- bad requirements;
- insufficient tests;
- unsound unsafe dependencies;
- supply-chain compromise;
- incorrect ground truth.

Therefore DNA uses a layered assurance model:

```text
type system / ownership
        ↓
compiler + lints
        ↓
unit + integration + property tests
        ↓
fuzz + mutation + adversarial validation
        ↓
real-trace scientific validation
        ↓
release provenance + operational evidence
```

Each layer checks a different class of invariant.

## Consequences

### Positive

- refactoring can fail at compile time when architectural assumptions are violated;
- coordinate, topology, strand, and state distinctions can become explicit domain vocabulary;
- malformed binary input is easier to contain behind a small validated parser boundary;
- concurrency can be added later without making data races an accepted risk in safe code;
- scientific evidence can remain immutable while derived interpretations evolve;
- the production binary can remain small and operationally simple;
- code review can focus more on scientific semantics because many structural mistakes become machine-detectable.

### Cost

- domain modeling requires more design before implementation;
- some changes require touching multiple exhaustive matches and constructors;
- compile times and borrow/type errors add friction during early experimentation;
- overusing newtypes, generics, or typestate could create unnecessary complexity;
- Python may remain faster for exploratory work and C++ may offer easier access to some established bioinformatics libraries.

The project should accept the useful friction while rejecting type-level complexity that does not eliminate a real failure mode.

## Design rule

When choosing between two implementations with comparable scientific behavior, prefer the design that makes the important invariant **explicit, local, typed, testable, and difficult to bypass**.

That is the primary reason DNA is implemented in Rust.
