---
implements: [SRS-NOM]
---

# Notation

When the [[reference/profiles|target profile]] recorded in `provenance.profile`
declares a notation style, the document carries a `notation` object; otherwise
the key is omitted. The examples below use the shipped human-mtDNA profile
(`human-mtdna-rcrs`), which pins the rCRS and declares `per_base_decimal`.

```json
"notation": {
  "style": "per_base_decimal",
  "calls": [
    {"call": "73G", "reads": ["A12_HV2F_03"]},
    {"call": "249DEL", "reads": ["A12_HV2F_03", "A12_HV3R_03"]},
    {"call": "309.1C", "reads": ["A12_HV3R_03"]},
    {"call": "315.1C", "reads": ["A12_HV3R_03"]}
  ]
}
```

## How it is produced

1. Each read's eligible variants (those with an empty exclusion list in
   `variants[].support[]`) are taken on their own.
2. They are normalized under the profile's indel placement (`right`: 3'/right-most,
   without crossing the reference origin seam).
3. The profile's nomenclature windows then name repeat haplotypes; the
   human-mtDNA profile names the HVS-II 303-315 and HVS-I 16180-16193 poly-C
   windows and the HVS-III 512-524 AC repeat in the lab's reviewer form
   (`309.1C 315.1C`, `310C 315DEL`, `513A 523DEL 524DEL`,
   `16183C 16184A 16189C`). A window that an edit straddles keeps its
   normalized form instead; the read's other windows are still represented.
4. The `per_base_decimal` style renders the represented variants one call per
   changed base:
   - substitution: `<position><base>`, for example `73G`;
   - deletion: `<position>DEL` for each deleted base, for example `523DEL 524DEL`;
   - insertion: `<anchor>.<ordinal><base>` after the preceding reference base,
     for example `309.1C 309.2C`; an insertion before position 1 uses anchor `0`.

5. A [[reference/profiles#schema-version-1|report convention]] of the profile labels
   each call it covers: `convention` gives its `name` and `label`
   (`not_recorded`: the reviewers do not record the call). The call is kept,
   and the haplotype is unchanged (SRS-NOM-017). A call that one read
   containing it does not label stays unlabelled. The human-mtDNA profile
   labels the length calls of the uninterrupted HVS-I C run that `16189C`
   creates, for example:
   `{"call": "16193.1C", "convention": {"name": "hvs1_c_length", "label": "not_recorded"}, "reads": ["HV1"]}`.

`calls` is ordered by position and then insertion ordinal. Each call lists, in
`reads[]` order, every read whose represented variants contain it.

## Interpretation boundary

`notation` makes sequence-equivalent descriptions from different reads converge
on one name; it is not a consensus. A read that covers a position without
reporting a call is not listed, and reads that disagree contribute different
calls. Use `coverage` and `locus_differences` to judge support and
disagreement. Ineligible observations never appear in `notation`.

Sanger repeat-artifact interpretation and primer callable ranges are not
applied.
