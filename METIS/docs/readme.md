```{image} metis_scopesim_logo.png
:width: 600px
:alt: METIS + ScopeSim
:align: center
```

# METIS + ScopeSim

## Introduction

The METIS data simulator is based on the generic simulator ScopeSim, a
descendant of the older SimCado/SimMETIS interface. METIS itself is handled
as an instrument package containing configuration files for the various
instrument modes and data files describing the instrument components.

The simulator currently supports imaging and long-slit spectroscopy modes.
The LM-band high-resolution IFU (LMS) mode is under development.

```{note}
**Bug reports and help desk**

If you come across a bug or get stuck with ScopeSim or the METIS package,
please [open an issue on GitHub](https://github.com/AstarVienna/irdb/issues)
or contact us by email (see below).

**Your feedback is the only way we know** what needs to be
changed or improved with the package and the simulator.

Please always include the output of `scopesim.bug_report()` from your
installation.
```

---

## Downloading the METIS instrument package

Once ScopeSim is installed, download the METIS instrument package into your
working directory:

```python
import scopesim
scopesim.download_packages(["Armazones", "ELT", "METIS"])
```

This installs the packages into the subdirectory `./inst_pkgs/`.

Your working directory should look like this afterwards:

```
my_simulations/
├── <your notebook>.ipynb
└── inst_pkgs/
    ├── Armazones/
    ├── ELT/
    └── METIS/
        └── docs/
            └── example_notebooks/
                └── <notebook>.ipynb   ← copy to working dir before running
```

---

```{include} ../../docs/ScopeSim_guide.md
```

---

## Example notebooks

Find the notebooks locally in `inst_pkgs/METIS/docs/example_notebooks/`
after downloading the package, or download them from the
[GitHub repository](https://github.com/AstarVienna/irdb/tree/dev_master/METIS/docs/example_notebooks).

```{warning}
Run notebooks in your working directory, **not** inside `inst_pkgs/`.
Copy the desired notebook out first.
```

### Introductory notebooks

| Notebook | Description |
|---|---|
| [`Introduction_to_Scopesim_for_METIS.ipynb`](example_notebooks/Introduction_to_Scopesim_for_METIS.ipynb) | Introductory overview of how to run simulations in ScopeSim. Also available as [PDF](example_notebooks/Introduction_to_Scopesim_for_METIS.pdf). |

### Scientific use-case notebooks

| Notebook | Description |
|---|---|
| [`IMG_L_N-examples.ipynb`](example_notebooks/IMG_L_N-examples.ipynb) | Imaging of HL Tau and an AGN model in the L and N bands |
| [`LSS-YSO_model_simulation.ipynb`](example_notebooks/LSS-YSO_model_simulation.ipynb) | Long-slit spectroscopy in the L-band of a young stellar object model |
| [`LSS_AGN-01_preparation.ipynb`](example_notebooks/LSS_AGN-01_preparation.ipynb) + [`LSS_AGN-02_simulation.ipynb`](example_notebooks/LSS_AGN-02_simulation.ipynb) | Long-slit spectroscopy in the N-band of an AGN model. Part 1: input data preparation. Part 2: simulation. |

### Effect demonstration notebooks

These notebooks are in `docs/example_notebooks/demos/`.

| Notebook | Description |
|---|---|
| [`demo_adc_wheel.ipynb`](example_notebooks/demos/demo_adc_wheel.ipynb) | Using the atmospheric dispersion correctors |
| [`demo_auto_exposure.ipynb`](example_notebooks/demos/demo_auto_exposure.ipynb) | Selecting `dit`/`ndit` automatically |
| [`demo_chopping_and_nodding.ipynb`](example_notebooks/demos/demo_chopping_and_nodding.ipynb) | Chop-nod difference images in the N band |
| [`demo_detector_modes.ipynb`](example_notebooks/demos/demo_detector_modes.ipynb) | Setting detector readout modes |
| [`demo_filter_wheel.ipynb`](example_notebooks/demos/demo_filter_wheel.ipynb) | Using the filter wheel(s) |
| [`demo_lss_simple.ipynb`](example_notebooks/demos/demo_lss_simple.ipynb) | Basic long-slit spectroscopy procedure |
| [`demo_grating_efficiency.ipynb`](example_notebooks/demos/demo_grating_efficiency.ipynb) | Spectral (grating) efficiency |
| [`demo_slit_wheel.ipynb`](example_notebooks/demos/demo_slit_wheel.ipynb) | Using the slit wheel for spectroscopy and imaging |
| [`demo_rectify_traces.ipynb`](example_notebooks/demos/demo_rectify_traces.ipynb) | Wavelength-calibrated and rectified 2D spectra |

### Instrument homepage

[METIS at Leiden Observatory](https://metis.strw.leidenuniv.nl/)
