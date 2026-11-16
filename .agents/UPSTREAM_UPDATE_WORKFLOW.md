# Upstream update workflow

Use this workflow whenever the fork needs to catch up with `upstream/dev`. The goal is to understand incoming changes, preserve the scientific intent of the custom overlays, and make only the compatibility updates that are needed.

## 1. Establish the starting point

- Read [the patch overview](PATCHES.md), the relevant overlay READMEs, and the open items in [TODO.md](TODO.md).
- Check the current branch, remotes, and working tree with `git status --short`, `git branch -vv`, and `git remote -v`.
- Preserve any existing user changes. Do not reset, overwrite, or stash them without direction.
- If the working tree or branch state makes a safe comparison unclear, report that before proceeding.

## 2. Fetch and identify upstream changes

Fetch updates the remote-tracking ref; it does not merge upstream into the fork:

```sh
git fetch upstream
```

Then establish the base and enumerate upstream commits and changed paths:

```sh
git merge-base HEAD upstream/dev
git log --oneline --decorate HEAD..upstream/dev
git diff --name-status "$(git merge-base HEAD upstream/dev)" upstream/dev
```

Inspect the relevant upstream commit diffs, not just commit subjects or file names. If the local remote-tracking ref was already current, say so and continue using that ref. Do not infer the contents of new upstream commits from an outdated local ref.

## 3. Trace interactions with the overlays

For every changed upstream file:

1. Check whether a corresponding copy or related implementation exists under `patch/tepper` or is selected by a custom Makefile.
2. Compare the incoming upstream version with both the fork's main-tree version and its overlay copy.
3. Read the relevant overlay README and follow cross-package dependencies documented in [PATCHES.md](PATCHES.md).
4. Separate upstream behavior from fork-specific behavior. In particular, do not assume a whole overlay file is custom; some overlays are synchronized source copies.
5. Consider interactions among combined variants such as Oagertz + Nexus/DICE + restart + AGAMA, not only each file in isolation.

Record which overlays are unaffected, which need synchronization, and which have a possible behavior conflict. Call out unresolved provenance or scientific uncertainty instead of guessing.

## 4. Report findings and propose changes

Before changing code, summarize:

- the upstream commits and relevant files reviewed;
- which overlay groups and build recipes they touch;
- compatibility concerns, conflicts, and any open questions;
- the smallest proposed changes and how each preserves the custom scientific behavior.

Follow the project instructions in the root [AGENTS.md](../AGENTS.md): do not modify code, merge, or rebase without approval. If the request already authorizes a specific implementation, work within that scope; otherwise, wait for approval before making those changes.

## 5. Update overlays and build integration

After the implementation is authorized:

- Update only the necessary overlay files, carefully carrying forward both upstream changes and the intended custom behavior. Avoid wholesale replacement when it would discard custom logic.
- Update affected custom Makefiles for source additions, module ordering, compiler options, and selected patch combinations.
- Keep unrelated refactoring out of the upstream sync. Track cleanup work separately in `TODO.md`.
- Update the affected patch README and [PATCHES.md](PATCHES.md) when behavior, dependencies, or known issues change. Preserve the README audit work as a distinct maintenance task where broader review is needed.

## 6. Integrate upstream

Once compatibility changes have been reviewed, choose and approve an integration strategy. This fork has a history of merging `upstream/dev`; merging is generally the straightforward way to preserve that history. Do not rebase the published `dev` branch unless explicitly directed. Resolve conflicts with the overlay analysis in view, and inspect the resulting tree and history before considering integration complete.

## 7. Verify when authorized

Do not run tests, builds, or simulations automatically. Agree on suitable checks for the affected patch combinations and representative models, then run them only when authorized. Report exactly what was run, the result, and any checks not performed. Compilation alone does not establish scientific correctness; state remaining uncertainty plainly.

## 8. Close out

Summarize the upstream range integrated, overlays and Makefiles changed, compatibility decisions, documentation updates, verification results, and remaining issues. Do not commit or push without explicit approval. If a commit was created but not pushed, remind the user that it remains local.
