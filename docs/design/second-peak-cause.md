# Second-peak cause

A call's second peak (SRS-CALL-015) has several causes. A reviewer tells them
apart by the peak's shape against the call's neighbours. `dna-sanger`'s
`second_cause` module labels each call with one of them (SRS-CALL-021). No
verdict, mask, consensus rule, or review item reads the label. It explains a
second peak to a reviewer, and the evidence adapter turns it into three
generic assertions the core may read (SRS-CALL-022, ADR-0080):
- whether the second base stands out;
- whether the call shows a second peak it does not explain;
- whether an isolated mixed call whose peak is a neighbour's tail or a dye
  blob may confirm the reference.

## Inputs

The rule reads the read's own calls only, never a reference position, an
alignment, or the vendor's peak positions:

- The second base and its share are the call's strongest other channel and
  its share of the two strongest peaks. This is the same computation that
  gives `CallEvidence.second`, so a report's share and label always agree.
- The call's centre is its primary channel's selected peak. Its spacing is
  half the distance between its neighbours' primary peaks.
- A channel's width is the width signal processing measures at the call's
  primary event: the samples around it where the channel stays at or above
  half its corrected amplitude. For the second channel this is its width
  around the primary event, not around its own apex. That is close for a
  centred peak, which is the only kind the width rule judges.

## Rule

A call is labelled when its second share exceeds
`signal_processing.second_floor`, the one second-channel floor of PROP-0009.
The first rule that holds gives the label:

| Label | Shape |
|---|---|
| `noise` | below `second_cause_noise_ratio` (3.0) times the call's lowest channel |
| `neighbour_tail` | no local maximum in the call's window, or a peak at least `second_cause_tail_offset` (0.4) spacings off the centre towards a neighbouring call of the second base |
| `off_centre` | that far off the centre, towards no call of its base |
| `dye_blob` | at least `second_cause_blob_width` (1.8) times as wide as the median primary peak of the `second_cause_width_calls` (6) calls on each side |
| `phase_shift` | centred, of a neighbouring call's base |
| `second_template` | centred, of another base, unexplained |

`run_signal::refine` turns a `second_template` into `context` when a run
signal holds the call's key at a usable run level and the share does not
exceed the call's threshold (SRS-CALL-020): a shadow the run's reads share.
Without a run signal, the shape alone decides.

A call whose callability state is mixed but whose own second share does not
exceed the floor is a `mixed_stretch`: its segment is mixed, not the call.

## Evidence

The research scan behind the rule covered raw traces of every current read,
over the whole reportable regions: 1.25 M calls with a second share of at
least 0.15. Among the calls DNA withholds today as mixed, the share whose
primary base would be wrong if they were lifted was:

| Label | Wrong |
|---|---|
| `context` | 0 % |
| `neighbour_tail` | 0.5 % |
| `dye_blob` | 1 % |
| `phase_shift` | 13 % |
| `second_template` | 24 % |
| `noise` | 27 % |
| `mixed_stretch` | 14 % |

Truth heteroplasmy is 0–0.04 % in every class, and most truth heteroplasmies
show a second share below 0.15. So the label is not a heteroplasmy detector;
heteroplasmy stays the review's two-strand rule.

Lifting the mixed verdict for `neighbour_tail` and `dye_blob` outright failed
the review gate on the development panel. The reference confirmation of
ADR-0080 is the narrow form that holds it.

`context` lifts nothing. It needs a run signal, which the consensus does not
read.
