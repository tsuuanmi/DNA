# Dependency Policy

Dependencies and lockfiles are part of the production codebase.

Before adding or materially changing a dependency:

1. confirm the capability is not already available;
2. justify maintenance, security, licensing, and supply-chain cost;
3. use the repository's package manager and lockfile workflow;
4. review the resulting lockfile diff;
5. update security/release evidence when the dependency changes the trust boundary.

Production dependencies must not be added speculatively. Release dependency review is tracked by [supply-chain security](../security/supply-chain.md).
