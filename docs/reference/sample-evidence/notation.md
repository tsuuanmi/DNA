# Human-mtDNA notation

When the reference is the rCRS (`provenance.reference.sha256` equals the
normalized rCRS sequence SHA-256
`f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66`), the
document carries a `notation` object. Against any other reference the key is
omitted.

```json
"notation": {
  "policy": "rcrs_right_aligned_hvs2",
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
2. They are right-aligned under the human-mtDNA policy, without crossing the
   rCRS origin seam.
3. The HVS-II 309/315 poly-C rule then expresses movement of the T310 anchor as
   C-run length changes. A read with an edit that straddles the validated
   303-315 window keeps its right-aligned form instead.
4. The represented variants are rendered one call per changed base:
   - substitution: `<position><base>`, for example `73G`;
   - deletion: `<position>DEL` for each deleted base, for example `523DEL 524DEL`;
   - insertion: `<anchor>.<ordinal><base>` after the preceding reference base,
     for example `309.1C 309.2C`; an insertion before position 1 uses anchor `0`.

`calls` is ordered by position and then insertion ordinal. Each call lists, in
`reads[]` order, every read whose represented variants contain it.

## Interpretation boundary

`notation` makes sequence-equivalent descriptions from different reads converge
on one name; it is not a consensus. A read that covers a position without
reporting a call is not listed, and reads that disagree contribute different
calls. Use `coverage` and `locus_differences` to judge support and
disagreement. Ineligible observations never appear in `notation`.

The HVS-III 513-524 and HVS-I 16189/16193 representation families, Sanger
repeat-artifact interpretation, and primer callable ranges are not applied.
