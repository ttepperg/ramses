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

Inspect the relevant upstream commit diffs, not just commit subjects or file names. Keep this changed-path list as the authoritative inventory for overlay discovery; retain it in the update notes before merging, since the comparison range changes after the merge. If the local remote-tracking ref was already current, say so and continue using that ref. Do not infer the contents of new upstream commits from an outdated local ref.

## 3. Trace interactions with the overlays

For every path in the changed-path inventory, search all of `patch/tepper` for corresponding overlay files, regardless of whether a current Makefile references them. Match the relative source path where possible; when overlays use different directory layouts, search by basename and relevant routine or symbol names, then verify each candidate by content and purpose. Do not use the current Makefiles to filter this discovery step.

For each candidate overlay:

1. Compare the incoming upstream version with the fork's main-tree version and the overlay copy.
2. Read its README and follow cross-package dependencies documented in [PATCHES.md](PATCHES.md).
3. Separate upstream behavior from fork-specific behavior. In particular, do not assume a whole overlay file is custom; some overlays are synchronized source copies.
4. Consider interactions among combined variants such as Oagertz + Nexus/DICE + restart + AGAMA, not only each file in isolation.
5. Record whether the candidate is active in a current Makefile, appears only in a legacy or alternate overlay, or has no current Makefile reference. Use Makefile `PATCH` order after discovery to establish which copy each configuration selects and to design build coverage; an absent or outdated Makefile must not hide an overlay candidate.

For each affected source, record the upstream change to carry forward, all corresponding overlay candidates found, the custom routines/data/conditional structure to preserve, and which Makefiles select each candidate. Treat patch READMEs and other documented references as evidence of custom intent; do not replace an overlay wholesale when it contains research-specific behavior. Preserve any user-identified structures pending later review, including the split `if (star.or.sink)` and `if (metal)` conditionals.

Ignore the contents of `patch/tepper/oagertz/src_orig` and `patch/tepper/namelist`; they are outside the tracked overlay review.

Record which overlays are unaffected, which need synchronization, and which have a possible behavior conflict. Call out unresolved provenance or scientific uncertainty instead of guessing.

## 4. Report findings and propose changes

Before changing code, summarize:

- the upstream commits and relevant files reviewed;
- which overlay groups and build recipes they touch;
- compatibility concerns, conflicts, and any open questions;
- the smallest proposed changes and how each preserves the custom scientific behavior.

Follow the project instructions in the root [AGENTS.md](../AGENTS.md): do not modify code, merge, or rebase without approval. If the request already authorizes a specific implementation, work within that scope; otherwise, wait for approval before making those changes.

## 5. Integrate upstream (first commit)

After the upstream changes and compatibility concerns have been reviewed and the integration is authorized, merge `upstream/dev` into `dev` as a distinct commit. Keep this commit focused on synchronizing the fork's main source tree with upstream; do not include patch-overlay synchronization in it.

- Preserve existing fork-specific changes while resolving conflicts; report uncertain behavior rather than silently discarding it.
- Review the merge result and commit it separately from all overlay changes.

Do not rebase the published `dev` branch unless explicitly directed. This fork has a history of merging `upstream/dev`, so a merge is generally the straightforward way to preserve that history.

## 6. Update overlays and build integration (second commit)

After the upstream merge, update overlays against the newly synchronized source. Keep this work in a separate commit from the upstream merge:

- Update only the necessary overlay files, carefully carrying forward both upstream changes and the intended custom behavior. Avoid wholesale replacement when it would discard custom logic.
- For the DICE mirror in `patch/tepper/dice`, review [`src_files.txt`](../patch/tepper/dice/src_files.txt), then preview and run [`sync.sh`](../patch/tepper/dice/sync.sh) from any directory with `bash patch/tepper/dice/sync.sh --dry-run` and `bash patch/tepper/dice/sync.sh`. The script resolves the repository root, checks all listed sources before copying, and uses the manifest paths as source data rather than shell commands.
- Update affected custom Makefiles for source additions, module ordering, compiler options, and selected patch combinations.
- Keep unrelated refactoring out of the upstream sync. Track cleanup work separately in `TODO.md`.
- Update the affected patch README and [PATCHES.md](PATCHES.md) when behavior, dependencies, or known issues change. Preserve the README audit work as a distinct maintenance task where broader review is needed.
- Review and commit the overlay changes separately from the upstream merge.

## 7. Verify when authorized

Once both commits are in place, do not run tests, builds, or simulations automatically. Agree on suitable checks for the affected patch combinations and representative models, then run them only when authorized. Report exactly what was run, the result, and any checks not performed. Compilation alone does not establish scientific correctness; state remaining uncertainty plainly.

Choose a small coverage matrix from the discovered overlay candidates and the actual `PATCH` combinations in the custom Makefiles. Include a broad integration configuration and separate configurations for affected candidates that the broad one does not select. Makefiles share object and module outputs in `bin`; clean between configurations or use isolated build trees so stale artifacts do not hide configuration-specific errors. For each smoke run, record the Makefile, model, duration/steps, and relevant feature flags. A passing run only covers optional code paths that were enabled and reached.

When verification runs on another system, record the commands and outcomes as user-reported results; do not describe them as locally run. Keep build/smoke evidence distinct from scientific validation and leave unresolved paths open in `TODO.md`.

## 8. Close out

Summarize the upstream range integrated, overlays and Makefiles changed, compatibility decisions, documentation updates, verification results, and remaining issues. Do not commit or push without explicit approval. If a commit was created but not pushed, remind the user that it remains local.
