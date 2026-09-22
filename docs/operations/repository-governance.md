# Repository Governance

This document defines the GitHub repository settings required for a production-ready
DNA release. These controls complement CI; they prevent required checks from being
bypassed.

## Default-branch ruleset

Protect `main` with a GitHub ruleset that applies to administrators as well as
contributors.

Required settings:

- require changes through pull requests;
- require at least one approving review;
- require review from Code Owners;
- dismiss stale approvals when new commits are pushed;
- require all review conversations to be resolved;
- require the branch to be up to date before merge;
- require linear history;
- block force pushes;
- block branch deletion;
- do not allow routine bypass of the ruleset.

Required pull-request status checks:

- `CI success` — aggregate gate for workflow security, Rust quality, dependency policy/review, MSRV, and repository/Python policy;
- `CodeQL Rust`;
- `ABIF fuzz smoke`.

The aggregate job must remain dependent on every mandatory job in `.github/workflows/ci.yml`; adding a new mandatory CI job requires adding it to `CI success` in the same change.

`OpenSSF Scorecard` runs on `main` and on a schedule and is therefore an
observability/security-posture signal rather than a pull-request merge check.

## Release-tag ruleset

Protect tags matching `v*`.

- only maintainers should be allowed to create or update release tags;
- release tags must not be force-updated or deleted as a normal operation;
- a tag must point to a commit reachable from `main`;
- the tag name must exactly match the crate version.

The release workflow independently verifies the final two conditions.

## Security settings

Enable:

- GitHub private vulnerability reporting;
- Dependabot alerts;
- Dependabot security updates;
- dependency graph;
- secret scanning and push protection when available for the repository plan.

The tracked `SECURITY.md` is the disclosure contract. Security-sensitive reports
must not be redirected into public issues.

## Merge strategy

Prefer squash merge for ordinary pull requests so `main` has one reviewed,
auditable change per PR. Disable merge commits if the repository ruleset enforces
linear history.

## Ownership

`.github/CODEOWNERS` declares repository-wide ownership and explicitly covers
high-risk release, dependency, contract, ADR, and production-source surfaces.

## Why settings are separate from workflows

A workflow can describe a required check but cannot make itself mandatory. The
ruleset is the enforcement layer that turns CI from an advisory signal into a
merge invariant.


## Bootstrap checklist

Repository files cannot enable GitHub security settings or rulesets by themselves. Before calling a release production-ready, an administrator must complete these GitHub settings:

1. **Settings → Security and quality → Advanced Security**
   - enable Dependency Graph;
   - enable Dependabot alerts;
   - enable Dependabot security updates;
   - enable Private vulnerability reporting;
   - enable Secret scanning and Push protection when available for the repository plan.
2. **Settings → Rules → Rulesets**
   - create the protected-`main` ruleset described above;
   - create the protected-`v*` tag ruleset described above;
   - keep bypass permissions empty for normal development.
3. **Settings → General → Releases**
   - enable release immutability so published release tags and assets cannot be replaced.
4. Re-run the pull request checks after Dependency Graph is enabled. The `Dependency review` job is intentionally required and will fail while the graph is disabled.

These settings are part of the production contract, not optional repository polish.
