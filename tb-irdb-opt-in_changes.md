# tb-irdb-opt-in

Four lines of YAML. They must merge together with the matching ScopeSim
branches, or ZShooter and X-shooter silently produce different numbers.

## Why

Two ScopeSim changes move behaviour that was previously unconditional behind a
setting that defaults to upstream:

- `!SIM.spectral.integrate_spectral_bins` (ScopeSim branch
  `rework/tb-05-integrate-spectral-bins`) — flux-conserving bin integration of
  spectra instead of upstream's point sampling.
- `fit_normalize` on the trace effects (ScopeSim branch
  `rework/tb-04-transform2d-normalize`) — rescaling `Transform2D` fit inputs to
  `[-1, 1]`.

Both packages were relying on those being the unconditional behaviour. Without
this branch they revert to upstream sampling and fits.

## Changes

| File | Line |
|---|---|
| `ZShooter_v2/default.yaml` | `integrate_spectral_bins: true` in the `alias: SIM` block |
| `XSHOOTER/default.yaml` | same |
| `ZShooter_v2/ZShooter_SPEC.yaml` | `fit_normalize: true` on `trace_list_analytical` |
| `XSHOOTER/XSHOOTER_SPEC.yaml` | `fit_normalize: true` on `spectral_traces` |

Both additions follow patterns the packages already use. Each `SIM` block
already carries `trace_flux_jacobian: true` for exactly this reason, and
`ZShooter_SPEC.yaml` already overrides `fit_degree: 6` on the same trace effect
that now gains `fit_normalize`.

## Merge-order hazard

This is the important part. Measured on a ZShooter `SPEC` observation, one
point source, fixed seed, 0.55-0.60 um:

| Configuration | plane 1 sum | plane 1 peak | SHA-256 |
|---|---:|---:|---|
| ScopeSim `main` + irdb `dev_master` | 1.117810243110e+04 | 4.4306 | `1117340b6f3ee3e8` |
| ScopeSim reworked + **this branch** | 1.117810243110e+04 | 4.4306 | `1117340b6f3ee3e8` |
| ScopeSim reworked + irdb `dev_master` | 1.152511427729e+04 | 5.5144 | `010a897c70bc3f0a` |

With this branch the simulation is **bit-identical** across all six image
planes. Without it, total flux moves +3.1% and the peak +24%, with **no error
and no warning** — the simulation simply runs and gives different answers.

So: merge this at the same time as, or before, the ScopeSim branches. Never
after.

`fit_normalize` is the less critical of the two. Measured separately, the
rescaling is immaterial for the centred focal-plane mm these traces use (both
routes reach machine precision), so it is set for exactness rather than
because it changes anything detectable. `integrate_spectral_bins` is the one
that carries the flux difference above.

## No other package is affected

`integrate_spectral_bins` and `fit_normalize` both default to the upstream
behaviour in ScopeSim, so METIS, MICADO, MOSAIC and the rest are unchanged by
either ScopeSim branch and need no edit here.

## Not in scope

`MANIFEST.in` is dirty in this working tree from unrelated earlier work and is
deliberately not part of this branch.
