# ZShooter detector engineering configurations

`SPEC` selects tapered qCCDs in B/G/R and 15 um LmAPDs in YJ/H/K.
`SPEC_ALTDETECT` selects untapered CCDs in B/G/R, the same YJ LmAPD, and 18 um H2RGs in H/K. Detector dimensions remain 2048 x 4096 in B/G/R and 2048 x 2048 in YJ/H/K. Optical prescription, focal-plane offsets and selected orders are shared.

Select one mode when constructing `UserCommands`, then construct an `OpticalTrain`. Reconstruct the train when changing the detector package or its initialization parameters, including pitch and QE. Exposure parameters can be supplied independently per arm through `readout(dit_blue=..., ndit_blue=..., ...)`.

The `!DET` operating points use input-referred electrons. LmAPDs have avalanche gain 10, usable well 15000 e, effective read noise 0.5 e RMS and photon-variance multiplier F=1.2. Their physical QE uses `H2RG_TER.dat`, the existing provisional NIR QE estimate; `APD_TER.dat` is the historical QE/F proxy and is no longer active. The adopted dark term remains independent of F. Changing F does not change the configured usable well or divide the already input-referred read noise again.

Other usable wells default to 60000 e. Read noise is 0.1 e for qCCDs, 2.5 e for CCDs and 6 e for H2RGs per independent integration. Dark rates are 0.5976, 2, 0.36 and 23 e/hour/pixel for qCCD, CCD, LmAPD and H2RG respectively. These are engineering specifications. The settings in `ZShooter_DET.yaml` and `ZShooter_ALTDETECT.yaml` are individually editable per detector.

DIT/NDIT are observing settings: no effect changes them automatically. `ExposureLimit` warns if the expected source+sky+thermal+dark charge per integration exceeds the usable well, after CCD charge summation and before random noise or software binning. It records `PEAKDIT` and `WELLDEP` without clipping. An integration means an independent exposure or completed ramp, not a nondestructive sample. The current CCD charge-packet limit is conservatively the same 60000 e as the native-pixel limit.

Native headers record `DETMODE`, `DETTYPE`, `PIXEL_SIZE` (mm), `EXNOISE` (variance factor), `AVGAIN` and `WELLDEP`, alongside actual exposure and detector WCS. Image charge stays input-referred; avalanche gain is distinct from the electronic electron/ADU conversion.

Analytical traces resolve their pixel pitch from the same detector commands. Seven along-slit samples and a sixth-degree transform fit keep the common physical footprint consistent across the two detector sizes. The larger H/K area changes coverage and sampling without retuning optics or reinstating the omitted YJ order.
