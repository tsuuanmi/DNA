# Consensus model

Owns `consensus`: consensus segments, decided sites with their run evidence,
and decision counts. `RunEvidence` answers the three questions its users ask
of a run's length evidence: whether the length is an estimate
(`is_estimate`), whether a read shows both of its ends (`shown_by_read`), and
whether the two-strand rule may decide it instead (`yields_to_two_strand`).

See [dna-consensus](../README.md).
