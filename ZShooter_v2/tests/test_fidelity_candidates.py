"""Physical checks using the configured ZShooter prescription and detector effects."""
from pathlib import Path
import numpy as np
import pytest
import yaml
from astropy import units as u
from astropy.io import fits
from types import SimpleNamespace
import scopesim as sim
from scopesim.effects import EchelleSpectralTraceList
from scopesim.effects.effects_utils import make_effect
from scopesim.effects.spectral_trace_list_utils import XiLamImage
from scopesim.detector import Detector
from scopesim.optics.image_plane_utils import header_from_list_of_xy

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def traces():
    config = yaml.safe_load((ROOT / 'ZShooter_SPEC.yaml').read_text())
    effect = next(item for item in config['effects'] if item['name'] == 'trace_list_analytical')
    kwargs = dict(effect['kwargs'])
    kwargs['filename'] = str(ROOT / kwargs['filename'])
    simulation = next(item for item in yaml.safe_load_all((ROOT / 'default.yaml').read_text()) if item['alias'] == 'SIM')
    kwargs['trace_flux_jacobian'] = simulation['properties']['spectral']['trace_flux_jacobian']
    return EchelleSpectralTraceList(**kwargs).spectral_traces


@pytest.mark.parametrize('channel', ['b', 'g', 'r', 'yj', 'h', 'k'])
def test_physical_slit_dispersion_across_sky_grids(traces, channel):
    trace = next(trace for name, trace in traces.items() if name.startswith(channel + '_'))
    assert 'optical_resolving_power' not in trace.meta
    wave0 = (trace.wave_min + trace.wave_max) / 2
    derivative = trace.slit_dispersion(wave0)
    assert derivative < 0
    sigma_wave, sigma_eta = wave0 / 100000, .65 / np.sqrt(8*np.log(2))
    wave = np.linspace(wave0 - .001*wave0, wave0 + .001*wave0, 2001)
    eta_ref = np.linspace(-.35, .35, 10001)
    illumination = np.exp(-.5*(eta_ref/sigma_eta)**2)
    eta_variance = np.trapezoid(eta_ref**2 * illumination, eta_ref) / np.trapezoid(illumination, eta_ref)
    expected_variance = sigma_wave**2 + derivative**2 * eta_variance
    for rows in [5, 10, 20]:
        step = .7 / rows
        eta = (np.arange(rows) + .5) * step - .35
        cube = np.exp(-.5*((wave-wave0)/sigma_wave)**2)[:, None, None] * np.exp(-.5*(eta/sigma_eta)**2)[None, :, None] * np.ones((1, 1, 3))
        header = fits.Header({'CTYPE1': 'LINEAR', 'CTYPE2': 'LINEAR', 'CTYPE3': 'WAVE', 'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec', 'CUNIT3': 'um', 'CDELT1': .1, 'CDELT2': step, 'CDELT3': wave[1]-wave[0], 'CRPIX1': 1., 'CRPIX2': 1., 'CRPIX3': 1., 'CRVAL1': 0., 'CRVAL2': eta[0], 'CRVAL3': wave[0]})
        fov = SimpleNamespace(cube=fits.ImageHDU(cube, header), meta={'xi_min': 0*u.arcsec, 'eta_min': -.35*u.arcsec, 'eta_max': .35*u.arcsec})
        image = XiLamImage(fov, lambda wave: 0., slit_dispersion=trace.slit_dispersion)
        flux = np.trapezoid(image.image[1], wave)
        variance = np.trapezoid((wave-wave0)**2 * image.image[1], wave) / flux
        # Midpoint quadrature across a seeing-limited 0.7 arcsec slit, five or more rows.
        assert variance == pytest.approx(expected_variance, rel=.03)
        assert flux == pytest.approx(np.sqrt(2*np.pi)*sigma_wave*np.trapezoid(illumination, eta_ref), rel=.015)


@pytest.mark.parametrize('detector_id,arm,dark,ron', [(0, 'blue', .000166, .1), (1, 'green', .000166, .1), (2, 'red', .000166, .1), (3, 'yj', .0001, .5), (4, 'h', .0001, .5), (5, 'k', .0001, .5)])
def test_dark_current_participates_in_poisson_statistics(detector_id, arm, dark, ron):
    config = yaml.safe_load((ROOT / 'ZShooter_DET.yaml').read_text())
    names = ['exposure_integration_selector', 'dark_current_selector', 'shot_noise', 'readout_noise_selector']
    effects = [item for item in config['effects'] if item['name'] in names]
    assert len(effects) == 4
    cmds = sim.UserCommands(properties={f'!OBS.dit_{arm}': 3600., f'!OBS.ndit_{arm}': 1, '!SIM.random.seed': 126})
    header = header_from_list_of_xy([-4.5, 4.5], [-4.5, 4.5], .015, 'D')
    detector = Detector(header, cmds=cmds, id=detector_id)
    for item in effects:
        make_effect(item, cmds).apply_to(detector)
    assert detector.data.mean() == pytest.approx(dark*3600, abs=.005)
    assert detector.data.var() == pytest.approx(dark*3600 + ron**2, rel=.02)


@pytest.mark.parametrize('detector_id,arm,bx,by', [(0, 'blue', 1, 2), (1, 'green', 2, 1), (2, 'red', 2, 2)])
def test_visible_ccd_charge_binning_is_independent_per_arm(detector_id, arm, bx, by):
    config = yaml.safe_load((ROOT / 'ZShooter_DET.yaml').read_text())
    configured = next(item for item in config['effects'] if item['name'] == 'ccd_charge_binning')
    properties = {f'!DET.bin{axis}_{channel}': 1 for channel in ('blue', 'green', 'red') for axis in ('x', 'y')}
    properties.update({f'!DET.binx_{arm}': bx, f'!DET.biny_{arm}': by})
    cmds = sim.UserCommands(properties=properties)
    effect = make_effect(configured, cmds)
    for identity in range(6):
        header = header_from_list_of_xy([-.3, .3], [-.3, .3], .015, 'D')
        detector = Detector(header, cmds=cmds, id=identity)
        detector.data[:] = 7.
        effect.apply_to(detector)
        expected_shape = (40//by, 40//bx) if identity == detector_id else (40, 40)
        assert detector.data.shape == expected_shape
        assert detector.data.sum() == pytest.approx(7*40*40)


@pytest.mark.parametrize('channel', ['b', 'g', 'r', 'yj', 'h', 'k'])
def test_native_trace_photon_normalization_is_independent_of_sky_grid(traces, channel):
    trace = next(trace for name, trace in traces.items() if name.startswith(channel + '_'))
    center = (trace.wave_min + trace.wave_max) / 2
    # Native detector geometry comes from the same IRDB prescription as the trace.
    config = yaml.safe_load((ROOT / 'ZShooter_DET.yaml').read_text())
    detector_effect = next(item for item in config['effects'] if item['class'] == 'DetectorList' and item['kwargs']['image_plane_id'] == trace.meta['image_plane_id'])
    detector_header = make_effect(detector_effect).image_plane_header
    for step in [.1, .05, .025]:
        wave = np.linspace(center*.998, center*1.002, 2001)
        xi = np.arange(-2, 2 + step/2, step)
        rows = round(.7 / step)
        eta_step = .7 / rows
        eta = (np.arange(rows) + .5) * eta_step - .35
        cube = np.exp(-.5*((wave-center)/(center/15000))**2)[:, None, None] * np.exp(-.5*(eta/.2)**2)[None, :, None] * np.exp(-.5*(xi/.35)**2)[None, None, :]
        header = fits.Header({'CTYPE1': 'LINEAR', 'CTYPE2': 'LINEAR', 'CTYPE3': 'WAVE', 'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec', 'CUNIT3': 'um', 'CDELT1': step, 'CDELT2': eta_step, 'CDELT3': wave[1]-wave[0], 'CRPIX1': 1., 'CRPIX2': 1., 'CRPIX3': 1., 'CRVAL1': xi[0], 'CRVAL2': eta[0], 'CRVAL3': wave[0]})
        metadata = {'wave_min': wave[0]*u.um, 'wave_max': wave[-1]*u.um, 'xi_min': xi[0]*u.arcsec, 'xi_max': xi[-1]*u.arcsec, 'eta_min': -.35*u.arcsec, 'eta_max': .35*u.arcsec}
        fov = SimpleNamespace(cube=fits.ImageHDU(cube, header), header=header, detector_header=detector_header, meta=metadata, trace_id=trace.trace_id)
        image = trace.map_spectra_to_focal_plane(fov)
        expected = cube.sum() * step * eta_step * (wave[1]-wave[0])
        assert image.data.sum() == pytest.approx(expected, rel=.01)
