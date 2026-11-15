# To do

1. [ ] Track the origin of the only direct source difference from upstream I found outside those additions is in pm/feedback.f90: it multiplies supernova energy by f_esn in two places.
2. [ ] Check whether custom sink overlays are compatible with upstream sink settings; if not, try make to make it compatible
3. [ ] Ignore oagertz/src_orig
4. [ ] Refactor patches, removing unnecessary duplicates, e.g. read_params.f90, by expanding the patch's root read_params.f90; modify Makefiles accordingly
5. [ ] Audit the patch READMEs against their corresponding upstream source files. For each overlay, verify that the README accurately describes the custom differences from upstream (separating them from synchronized upstream code), and update it where needed. Include relevant dependencies/interactions and any known issues, without duplicating details unnecessarily.
6. [ ] Resolve the sink.and.star + metal conflicting declarations
7. [ ] Address bugs reported in READMEs

## Future tasks

0. [ ] Merge nexus with dice, and remove unnecessary patches
1. [ ] Fix dice_restart patch + star formation
2. [ ] MAJOR: Migration of patches to Ramses 2
