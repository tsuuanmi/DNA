//! Trace signal engineering (ADR-0084, PROP-0012): numerical work on a read's
//! trace windows that a modality flags, from plain data. The crate reads no
//! instrument format, reference, alignment, or truth; the modality turns its
//! results into kernel assertions.

pub mod decompose;
