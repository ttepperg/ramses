# RAMSES fork patches

This document summarizes the custom patch set in `patch/tepper` and how it relates to the upstream RAMSES source. It is a map of the current organization, not a file-by-file diff. The patch files are source overlays that preserve RAMSES's package layout; many contain both custom changes and upstream code copied forward so the custom behavior continues to apply.

## Upstream and comparison point

- Fork remote: `origin` → `git@github.com:ttepperg/ramses.git`.
- Upstream remote: `upstream` → `https://github.com/ramses-organisation/ramses.git`.
- Compared branch: local `dev` (`7981357`, after merging upstream), against local `upstream/dev` (`386c7742`, 9 November 2026).
- Relevant common ancestor: `386c774258c5316418c5038c429e8dcb91e08832` (`Add more namelist checks for hydro parameters (#341)`). This commit is the current local `upstream/dev` tip and is included in the fork's merge commit.

The branch history includes repeated merges from upstream. In the direct source comparison against `upstream/dev`, the only difference outside the patch and build additions is `pm/feedback.f90`: two supernova energy expressions are multiplied by `f_esn`. The current evidence establishes the code difference, but not its origin or whether it was intended as a fork-specific change; see [TODO.md](TODO.md), item 1. Do not assume it is a custom behavior until that is traced.

## Patch groups

### Output locations, units, and initialization

The `amr`, `hydro`, and `poisson` overlays add `output_dir` support to paths used for outputs, restart files, and initialization. `amr/units.f90` adds unit conversions and adjusts scaling for gravity runs. These changes touch multiple packages because output and initialization work is spread across AMR, hydro, particle, and Poisson routines. A build using this behavior needs the relevant overlays together; a partial selection may leave some files using the default output location.

See [amr/README](../patch/tepper/amr/README), [hydro/README](../patch/tepper/hydro/README), and [poisson/README](../patch/tepper/poisson/README).

### Particle, sink, and restart behavior

The `pm` overlay extends sink particle initialization and evolution: it adds optional sink mass damping, support for adding sinks at restart, input scaling, and sink CSV output. It also updates particle output and routes relevant files through `output_dir`. The current sync carries forward upstream's per-level AGN feedback reservoir and optional gas/particle dynamical-friction support while retaining those custom sink settings. Compatibility with the custom restart path and runtime behavior still needs validation; see [TODO.md](TODO.md).

The standalone `dice_restart` overlay adds a `ramses_dice` restart loader and can optionally add DICE-format initial-condition particles to a restart. It depends on particle, restart, and DICE/Nexus code paths and must be combined through the appropriate Makefile. The related `oagertz/dice_restart` is a combined variant built on the Oagertz and Nexus/DICE overlays.

The restart READMEs report unresolved concerns: stellar ages may be wrong when restart time is restored (and setting time back to zero has its own consequences), and adding initial-condition particles at restart may create duplicate particle IDs. They also note particles at or beyond the box boundary may be discarded. See the linked READMEs for details and current status.

See [pm/README](../patch/tepper/pm/README), [dice_restart/README](../patch/tepper/dice_restart/README), and [oagertz/dice_restart/README](../patch/tepper/oagertz/dice_restart/README).

### DICE and legacy Nexus

`patch/tepper/dice` contains copies of the DICE implementation now present in upstream RAMSES. Its README says these files are retained to support the fork's patch/build arrangement; it is not, by itself, a distinct custom DICE implementation.

`patch/tepper/dice-legacy` is the older DICE-derived loader and related code retained as the basis of the Nexus framework. It was renamed and its parameters, file types, and module names were changed to avoid collisions with upstream DICE. It includes custom initial-condition loading and related particle, gas, and refinement routines. The `oagertz/dice` and `oagertz/agama_dice` variants compose this loader with the Oagertz and AGAMA changes.

The DICE/Nexus arrangement is transitional: the TODO list proposes merging Nexus with DICE and removing redundant patches. Avoid treating today's folder boundaries as a final architecture.

See [dice/README](../patch/tepper/dice/README), [dice-legacy/README](../patch/tepper/dice-legacy/README), and the corresponding READMEs under [oagertz](../patch/tepper/oagertz/README).

### Oagertz physics and variants

`patch/tepper/oagertz` maintains an overlay based on Oscar Agertz's RAMSES changes. The central features are stellar feedback and star formation, including metal yields and the extra star-particle properties `mpb` (birth mass) and `zp` (metallicity). Supporting changes span particle initialization, communication, output, AMR and hydro parameters, and feedback. Its files are maintained against upstream RAMSES, so the overlay includes synchronized upstream code alongside the Oagertz-specific behavior.

Subdirectories select combinations or variants:

- `oagertz/dice` and `oagertz/dice_restart` combine Oagertz behavior with Nexus/DICE loading and restart support.
- `oagertz/agama_dice` combines the AGAMA and DICE paths.
- `oagertz/wh2007` selects Woosley & Weaver (2007) Type II supernova O/Fe yields in place of the default yield functions.
- `oagertz/src_orig` holds reference copies used in the documented synchronization workflow; it is not an active patch directory to compile.

The standalone and DICE Oagertz particle-tree overlays retain their custom particle-property communication while following upstream's deferred linked-list move updates. The shared Oagertz `pm_commons` declares the move-tracking arrays, allocated by `patch/tepper/amr/init_amr.f90`. The Oagertz hydro-parameter reader follows upstream's split solver/refinement readers while retaining its custom feedback and star-formation parameters.

Some overlay interactions are delicate. The READMEs document a `star`/`sink` plus `metal` conditional issue and warn to preserve allocation, loop, and preprocessor logic when updating files. Known Oagertz and restart issues are recorded in the subpatch READMEs and tracked at a high level in [TODO.md](TODO.md).

See [oagertz/README](../patch/tepper/oagertz/README), [oagertz/dice/README](../patch/tepper/oagertz/dice/README), [oagertz/dice_restart/README](../patch/tepper/oagertz/dice_restart/README), and [oagertz/wh2007/README](../patch/tepper/oagertz/wh2007/README).

### AGAMA gravitational potentials

`patch/tepper/agama` exposes an installed AGAMA library to RAMSES. It supports two approaches: `agama/rho` supplies an analytic density (`rho_ana`) for RAMSES to solve for the potential and force; `agama/grav` evaluates and supplies AGAMA forces directly. The reader supports scheduled potential-file changes using a time-to-file table, with a fixed-file fallback. `agama/interface` demonstrates a custom wrapper around the standard AGAMA Fortran interface, including explicit potential-object deletion, without rebuilding AGAMA itself.

This patch requires AGAMA to be built on the system where RAMSES will run. Select one AGAMA implementation in the Makefile and configure the matching Poisson settings; the two approaches are not interchangeable at runtime. `oagertz/agama_dice` supplies a composition for the corresponding Oagertz/DICE use case.

See [agama/README](../patch/tepper/agama/README), [agama/interface/README](../patch/tepper/agama/interface/README), and [oagertz/agama_dice/README](../patch/tepper/oagertz/agama_dice/README).

## Build integration

The custom `bin/Makefile_*` files select compatible overlay combinations and external options, including DICE/Nexus, restart loading, Oagertz feedback and yields, and AGAMA density or force modes. They are the integration layer for the patch groups above; use the Makefile indicated by the relevant model's configuration and patch README. The upstream default Makefile and the custom recipes have been periodically harmonized, so future upstream build changes may require updating these recipes.

See [bin/README](../bin/README) for the listed configurations and example model pairings.

## Maintenance notes

- The patch files are maintained source overlays, not necessarily minimal diffs. A file's presence here does not mean every line is custom; consult its README and compare against upstream before attributing behavior.
- Several overlays duplicate shared files such as `read_params.f90`. The TODO list calls for reducing these duplicates and updating Makefiles; dependencies described above reflect the current layout and may change during that refactor.
- Restart/star-formation behavior, sink compatibility, and the `star`/`sink` plus `metal` conditionals have open questions or reported issues. See [TODO.md](TODO.md) and the specific patch READMEs for details. This overview does not claim those issues are resolved.
- The TODO list also tracks broader documentation review and a future migration to RAMSES 2.

This is an exploratory summary. No builds or tests were run while preparing it.
