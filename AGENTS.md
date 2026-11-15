
# RAMSES Fork — Project Instructions

## Purpose
This repository is a fork of the upstream RAMSES codebase,
containing custom patches (patch/tepper) for scientific research.

The objective is to keep these patches compatible with
upstream developments while preserving their intended
scientific behaviour.

## Patch documentation
Read .aeon/PATCHES.md before proposing upstream updates.
This file documents the purpose, implementation, and
dependencies of the custom modifications, in addition to critcal unresolved issued (also documented in each patches README).

## Upstream update workflow
1. Inspect changes introduced upstream.
2. Identify interactions with our custom patches.
3. Explain any conflicts or compatibility concerns.
4. Propose the smallest necessary modifications.
5. Preserve the scientific intent of all custom patches.
6. Update PATCHES.md if the implementation changes.

## Working rules
- Never merge, rebase, or modify code without approval.
- Do not automatically run tests or simulations.
- Avoid unrelated refactoring.
- Prefer minimal, understandable changes.
- Report any uncertainty about scientific correctness.
- Do not commit or push without explicit approval.
