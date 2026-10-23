# XSHOOTER

VLT UT3 / Paranal slit spectroscopy for ScopeSim, with simultaneous UVB, VIS, and NIR detector images. The package requires the accompanying `xshooter` ScopeSim changes. It is a development model: five-percent end-to-end accuracy is not established, particularly in the NIR, telluric bands, order edges, and wide optical slits.

The companion notebook and full model/accuracy notes are in the zs-scopesim workspace:

- `notebooks/01_Kilonovae_XSHOOTER.ipynb`
- `docs/xshooter.md`

## Configuration

```python
cmds = scopesim.UserCommands(use_instrument='XSHOOTER', set_modes=['SPEC'])
train = scopesim.OpticalTrain(cmds)
train.observe(source)
uvb, vis, nir = train.readout()
```

Set the properties before constructing the optical train:

| Property | Values |
|---|---|
| `!INST.uvb_curr_slit` | `'0.5', '0.8', '1.0', '1.3', '1.6', '5.0'` |
| `!INST.vis_curr_slit` | `'0.4', '0.7', '0.9', '1.2', '1.5', '5.0'` |
| `!INST.nir_curr_slit` | `'0.4', '0.6', '0.9', '1.2', '5.0', '0.6JH', '0.9JH'` |
| `!OBS.dit_uvb`, `!OBS.dit_vis`, `!OBS.dit_nir` | Seconds per integration; NIR noise table covers 1–1800 s |
| `!OBS.ndit_uvb`, `!OBS.ndit_vis`, `!OBS.ndit_nir` | Number of integrations |
| `!DET.binx_uvb`, `!DET.biny_uvb`; likewise `_vis` | Pairs `(1, 1)`, `(1, 2)`, `(2, 2)` in FITS column/row order |
| `!DET.readout_mode_uvb`, `!DET.readout_mode_vis` | `'slow'` (100 kHz high gain), `'fast'` (400 kHz low gain) |
| `!OBS.seeing`, `!OBS.airmass`, `!ATMO.pwv` | Zenith seeing at 500 nm in arcsec, airmass, and PWV in mm |
| `!ATMO.local_sky_file` | `false` for online SkyCalc, or a physical ESO SkyCalc FITS file matching the conditions |

NIR has no on-chip binning. Slits are 11 arcsec long. Native arrays are UVB `(3000, 2048)`, VIS `(4000, 2048)`, and NIR `(1100, 2048)` in NumPy y/x order. These are complete standard science readout regions; optical prescan/overscan is omitted. No optional windows or IFU are implemented.

The three returned HDULists each contain an empty primary HDU and a detector ImageHDU. Images are bias-subtracted floating ADU. Optical NDITs sum; NIR NDITs average. `DIT`, `NDIT`, `EXPTIME`, `COMBINE`, gain, slit, arm, and binning are recorded. Readout overheads, saturation/TLI recovery, persistence, nonlinear ramp processing, bad-pixel maps, and fringing realizations are outside this version. Use exposures below saturation.

`Paranal_spectroscopy.yaml` and `VLT_UT3_Cassegrain.yaml` are separate components. A telescope substitution must replace the latter, not append a second telescope. The optical response calibration is conditional on the adopted telescope reflectivity and collecting area.

## Calibrations

`prepare_calibrations.py` builds the runtime trace, response, QE, NIR background, and NIR noise files from pinned inputs. Runtime simulation does not compile or load the ESO kernel.

```bash
conda run -n zssim python XSHOOTER/prepare_calibrations.py
conda run -n zssim python XSHOOTER/calibration/fit_lsf.py
conda run -n zssim python -m pytest XSHOOTER/tests
```

The trace tables contain 12 UVB, 15 VIS, and 16 NIR orders. Each has wavelength in micrometres, along-slit position in arcseconds, physical raw-oriented focal coordinates in millimetres, and a cross-slit wavelength derivative in micrometres/arcsecond. `FITDEG` specifies the fitted polynomial degree; `SLITCOL` selects the derivative; `LSFRES` supplies optical Gaussian resolving power. Slit integration, optical blur, and detector sampling remain distinct calculations.

`calibration/fit_lsf.py` reproduces the arm-level optical blur fit against User Manual Table 13 and writes `slit_resolution.ecsv`. The default-slit estimates differ from the published values by 0.3%, 1.8%, and 2.6%; the widest tabulated UVB/VIS slits differ by 7.0%/5.7%. These are uniformly illuminated slit estimates, not a point-source resolution guarantee. The mapping tests independently check photon conservation, centroid, and line variance through the configured traces at two numerical sky-grid spacings.

`calibration/eso_model.c` contains the mathematical functions from the ESO X-shooter 3.8.12 physical-model kernel, plus a small C interface. Its source and GPL-2.0-or-later license are identified in the file and `calibration/ESO_COPYING`. The pinned optimized model parameters are from the following public ESO archive products:

| Arm | Physical model product |
|---|---|
| UVB | `M.XSHOOTER.2024-01-03T13:03:57.940` |
| VIS | `M.XSHOOTER.2024-01-03T13:03:59.573` |
| NIR | `M.XSHOOTER.2024-01-03T13:03:56.276` |

Internal ESO model coordinates are not raw detector pixels. Preparation undoes the UVB half-turn and NIR CPL quarter-turn/crop. In NIR, raw half-pixel-centered coordinates are `(2044 - model_y, model_x + 20)`. The raw image retains 2048 columns and 1100 rows. Physical focal coordinates are centered on each array. The polynomial fit is checked against the optical kernel on independent wavelength/slit samples; every order meets 0.03 pixel RMS and 0.1 pixel maximum forward/inverse residual. `calibration/trace_residuals.ecsv` records the actual results.

The optical response uses the public EG21 (`CPD-69 177`) standard-star nod pair below. Rates are measured in a 3.6 arcsec spatial aperture after nod subtraction. The response divides observed electron rate by reference photon rate, atmospheric transmission, two adopted aluminium mirror reflectivities, and detector QE. Continuum bins span 32 native dispersion pixels. The A nod supplies the response; B supplies the repeatability check. The normalization retains finite-aperture losses and any residual transparency error in that observation. This is not an independent absolute-throughput validation.

| Local input filename | ESO product ID |
|---|---|
| `standard-uvb.fits` | `XSHOO.2024-01-01T01:14:01.322` |
| `standard-b-uvb.fits` | `XSHOO.2024-01-01T01:18:24.675` |
| `standard-vis.fits` | `XSHOO.2024-01-01T01:14:06.522` |
| `standard-b-vis.fits` | `XSHOO.2024-01-01T01:18:29.906` |
| `standard-nir.fits` | `XSHOO.2024-01-01T01:14:09.698` |
| `standard-b-nir.fits` | `XSHOO.2024-01-01T01:18:32.604` |
| `response-uvb.fits` | `M.XSHOOTER.2024-01-02T18:39:02.690` |
| `response-vis.fits` | `M.XSHOOTER.2024-01-02T18:39:10.320` |
| `response-nir.fits` | `M.XSHOOTER.2024-01-02T18:39:01.123` |
| `atmos-ext-uvb.fits` | `M.XSHOOTER.2019-04-04T09:50:14.503` |

Fetch each at `https://dataportal.eso.org/dataPortal/file/<product ID>`. Raw downloads use UNIX compress (`.Z`); decompress with `gzip -dc` before giving them the `.fits` names above. The response products contain the pipeline's `REF` standard-star flux in erg/cm2/s/Angstrom; the `RESPONSE` column itself is not used as raw instrumental throughput.

`sky-default-raw.fits` is an ESO SkyCalc spectrum with observatory Paranal, airmass 1.3, PWV 2.5 mm, new Moon (`moon_sun_sep=0`), wavelength range 300–2500 nm, fixed spectral resolution 200000, and no LSF. Other sky parameters use the SkyCalc defaults. Its molecular transmission selects clear continuum intervals; its total transmission is used only for the exploratory NIR measurements. Optical normalization uses the archived extinction curve. The pipeline reference flux products and raw inputs are not required at runtime.

To repeat the standard calibration after downloading those inputs into an explicit directory:

```bash
conda run -n zssim python XSHOOTER/calibration/measure_standard.py --data-dir /path/to/inputs
conda run -n zssim python XSHOOTER/prepare_calibrations.py
```

`optical_response.fits` is the frozen optical calibration. `standard_measurements.fits` retains the binned A and B measurements. `standard_repeatability.ecsv` reports the B/A bias and RMS. The optical second-nod RMS is below 5% for 26/27 orders; UVB order 13 has about 16% RMS including low-response edges. VIS order 17 has only four continuum bins. Interpolation over absorption regions and endpoint continuation are explicit approximations.

NIR remains based on the ESO ETC instrument envelope, empirical correction, QE, and a scalar sinc-squared blaze centered on the manual's Table 14. The standard measurements show that this cannot currently support a general 5% claim. The original ETC ASCII curves are frozen in `calibration/eso_etc`; they can be retrieved through the [ESO calibration database](https://www.eso.org/observing/etc/bin/simu/calibdisp?INS.NAME=XSHOOTER+INS.MODE%3Dnone).

`nir_background.fits` uses Table 12's four warm-background measurements, including the JH-filter reduction. It is a smooth cross-dispersion approximation derived from the configured traces and slit widths, with endpoint continuation, not a measured 2D stray-light map. Sky, telescope emission, and detector dark current are separate. Temperature variation of this empirical component is not calibrated. The graph-digitized JH cutoff has a stated passband continuation; localized leaks and fringes are omitted.

## References

- [ESO documentation](https://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc.html)
- [User Manual v10 / P114](https://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc/VLT-MAN-ESO-14650-4942_P114v1.pdf), tables 6–14 and Figure 13
- [Vernet et al. 2011](https://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc/xshooter_Vernet2011.pdf)
- [Pipeline manual 3.8.12](https://ftp.eso.org/pub/dfs/pipelines/instruments/xshooter/xshoo-pipeline-manual-3.8.12.pdf)
- [ESO pipeline source kit 3.8.12](https://ftp.eso.org/pub/dfs/pipelines/instruments/xshooter/xshoo-kit-3.8.12-3.tar.gz)
- [ESO SkyCalc](https://www.eso.org/observing/etc/skycalc/)
