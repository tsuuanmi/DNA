# Trace Integrity Evidence

Part of the canonical [signal-processing method](../signal-processing.md).

DNA retains the validated PLOC series as the current method's event-anchor
authority, but optional vendor PBAS/PCON cardinality is treated separately.
A vendor series may be shorter or longer without creating or deleting DNA
calls. Public integrity evidence records the PLOC count and, when present, each
vendor-series count.

For adjacent PLOC values DNA records minimum, median, and maximum spacing in
trace-sample units. One-locus inputs have no spacing summary.

Analyzed DATA values are signed 16-bit samples. Values exactly equal to the
signed-16-bit extrema are counted as `clipped_channel_samples`. This is an exact
representation-boundary observation, not a general artifact classifier.

For each PLOC locus, let:

~~~text
event_signal = sum(A/C/G/T corrected amplitudes at the refined event sample)
~~~

Among loci with positive event signal, DNA records:

~~~text
maximum_to_median_event_signal_ratio =
    max(event_signal) / median(event_signal)
~~~

using the existing six-decimal metric rounding. This ratio is deliberately not
thresholded into a dye-blob/high-amplitude-artifact label. It is evidence for
validation and review until a classifier has its own specification and truth
data.

Trace-integrity observations never change calls, SNR windows, quality, trim,
alignment, or variants.
