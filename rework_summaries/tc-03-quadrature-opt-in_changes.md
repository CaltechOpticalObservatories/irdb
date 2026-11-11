# tc-03-quadrature-opt-in

Two lines. Must merge together with ScopeSim branch
`rework/tc-03-quadrature-order`, which changes the default.

## Why

ScopeSim's per-pixel quadrature order becomes
`!SIM.spectral.pixel_quadrature_order` with **default 1** — upstream's single
sample at the pixel centre. Both packages were relying on the fork's 4-point
integration being unconditional, so without this they silently revert.

| File | Line |
|---|---|
| `ZShooter_v2/default.yaml` | `pixel_quadrature_order: 4` in the `alias: SIM` block |
| `XSHOOTER/default.yaml` | same |

Both go beside the existing `trace_flux_jacobian: true`, which is already
there for exactly this reason.

## Measured

With this branch, a ZShooter `SPEC` observation is byte-identical to `main`
across all 216 fingerprinted quantities. Without it, 24 differ — image-plane
and readout sums move by up to +0.03% in the NIR arms. Small, but a real
integration-vs-point-sampling difference, and silent: no error is raised.

Merge this at the same time as the ScopeSim branch, never after.

## Not in scope

`MANIFEST.in` is dirty in this working tree from unrelated earlier work and is
deliberately not part of this branch.
