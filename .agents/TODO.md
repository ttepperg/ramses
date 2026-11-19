# To do

1. [ ] Track the origin of the only direct source difference from upstream I found outside those additions is in pm/feedback.f90: it multiplies supernova energy by f_esn in two places.
2. [ ] Validate the custom sink overlays after carrying forward upstream's per-level AGN feedback reservoir and optional gas/particle drag settings. In particular, check the interaction with the custom sink restart path and `init_sink.f90`, which was previously kept out of sync; runtime behavior still needs verification.
3. [ ] Review the split `if (star.or.sink)` / `if (metal)` conditionals in the Oagertz patches and determine whether the split structure is needed. Preserve the existing splits until then; consult `src_orig` if useful.
4. [ ] Refactor patches, removing unnecessary duplicates, e.g. read_params.f90, by expanding the patch's root read_params.f90; modify Makefiles accordingly
5. [ ] Audit the patch READMEs against their corresponding upstream source files. For each overlay, verify that the README accurately describes the custom differences from upstream (separating them from synchronized upstream code), and update it where needed. Include relevant dependencies/interactions and any known issues, without duplicating details unnecessarily.
6. [ ] Resolve the sink.and.star + metal conflicting declarations
7. [ ] Address additional bugs and uncertain behavior reported in READMEs:
   - [ ] `pm`: Restarting with a sink when the previous run had no sinks reportedly fails if the expected `sink_XXXXX.csv` file is absent. The README says attempted changes may have fixed this; verify.
   - [ ] `dice_restart`: Check why two snapshots may be written immediately after restart when only one is expected.
   - [ ] `dice_restart`: Determine whether resetting `mcons` and `econs` to zero on restart is a bug or intended behavior.
   - [ ] `oagertz`: Check whether the `S*unit_out` variables are defined when `SNdiagnostics = .false.`; the README says they may not be.
   - [ ] `oagertz`: Verify whether `kinetic_feedback` should always use `zp(...,1)` for `ZSN`.
   - [ ] `agama`: Reproduce the reported `gravana.f90` restart failure when no potential file is present (model `mmw_m1_a_pcg_hyr`, restart at `output_00227`) and determine whether the `clean_stop` is expected.
   - [ ] `agama/interface`: Investigate the small reported memory increase (less than 0.1 MB per `agama_initfromfile()` call) despite deleting the potential object.
8. [ ] Automate upstream-update inventory and compile coverage: take the upstream changed-path list, find candidate overlay copies throughout `patch/tepper` independently of Makefile references, then map candidates to representative custom Makefiles by `PATCH` composition. Build in isolated trees or clean between configurations because `bin` artifacts are shared, and capture per-configuration results. Keep model-specific smoke simulations manual until representative remote test cases and feature settings are defined.

## Future tasks

0. [ ] Merge nexus with dice, and remove unnecessary patches
1. [ ] Migrate the restart functionality from the old, no-longer-maintained `patch/init/restart` into `dev`, merge it with the DICE patch, and retain the additional features in `patch/tepper/dice_restart`. The restart patch works except when stars are present: reading the restart time correctly can produce negative stellar masses, while resetting time to zero affects stellar ages and their subsequent evolution. Investigate the related `Main step`/`nstep_coarse` issue and the README's note about missing patched `init_time.f90`.
2. [ ] MAJOR: Migration of patches to Ramses 2
