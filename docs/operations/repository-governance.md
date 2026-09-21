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

- `Rust quality`;
- `Dependency policy`;
- `Dependency review`;
- `MSRV 1.88`;
- `Repository policy and Python tooling`;
- `CodeQL Rust`;
- `ABIF fuzz smoke`.

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
