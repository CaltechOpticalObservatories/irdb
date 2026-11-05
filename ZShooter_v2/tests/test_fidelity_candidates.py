"""Physical checks using the configured ZShooter prescription and detector effects."""
from pathlib import Path
from copy import deepcopy
import numpy as np
import pytest
import yaml
from astropy import units as u
from astropy.io import fits
from types import SimpleNamespace
import scopesim as sim
from scopesim.effects import EchelleSpectralTraceList
from scopesim.effects.effects_utils import make_effect
from scopesim.utils import from_currsys
from scopesim.effects.spectral_trace_list_utils import XiLamImage
from scopesim.detector import Detector
from scopesim.optics.image_plane_utils import header_from_list_of_xy

ROOT = Path(__file__).resolve().parents[1]


def test_science_defaults_draw_source_sky_and_dark_per_pixel(monkeypatch):
    monkeypatch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(ROOT.parent))
    ignore = ['atmo_transmission', 'continuum_emission', 'airglow_and_interline_continuum', 'adc_residuals']
    cmds = sim.UserCommands(use_instrument='ZShooter_v2', properties={'!OBS.mjdobs': 60310.25, '!SIM.random.seed': 38291})
    for document in cmds.yaml_dicts:
        if 'effects' in document:
            document['effects'] = [effect for effect in document['effects'] if effect['name'] not in ignore]
    train = sim.OpticalTrain(cmds)
    assert train['shot_noise'].include
    assert train['readout_noise_selector'].include
    # Different source and sky levels in nine regions, including pure dark,
    # sky-only and source-only pixels. No stochastic effect is disabled.
    # 90,000 pixels per region constrain even the sparse dark-only Poisson tail.
    source = np.repeat([0., .003, 3.], 300)[None, :]
    sky = np.repeat([0., .03, 5.], 300)[:, None]
    header = header_from_list_of_xy([-6.75, 6.75], [-6.75, 6.75], .015, 'D')
    for aid, arm in enumerate(('blue', 'green', 'red', 'yj', 'h', 'k')):
        train.cmds[f'!OBS.dit_{arm}'], train.cmds[f'!OBS.ndit_{arm}'] = 300., 2
        dark = train['dark_current_selector'].get_effect(aid)._get_dark_level(aid)
        ron = from_currsys(train['readout_noise_selector'].get_effect(aid).meta['noise_std'], train.cmds)
        bias = from_currsys(train['bias_selector'].get_effect(aid).meta['bias'], train.cmds)
        factor = train.detector_managers[aid][0]._native_header['EXNOISE']
        expected = (source+sky+dark)*600
        frames = []
        for roid in (0, 1):
            train.cmds['!OBS.roid'] = roid
            detector = Detector(header, cmds=train.cmds, id=aid)
            detector.data[:] = source+sky
            for effect in train.optics_manager.detector_effects:
                effect.apply_to(detector)
            residual = (detector.data-bias-expected)/np.sqrt((factor*(source+sky)+dark)*600+2*ron**2)
            for yy in (0, 300, 600):
                for xx in (0, 300, 600):
                    region = residual[yy:yy+300, xx:xx+300]
                    assert abs(region.mean()) < .05
                    assert region.var() == pytest.approx(1., abs=.06)
            frames.append(detector.data.copy())
        assert not np.array_equal(*frames)


def test_reused_configuration_matches_fresh_slits_and_exposures(monkeypatch):
    from synphot import SourceSpectrum, Empirical1D
    from synphot.units import PHOTLAM

    monkeypatch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(ROOT.parent))
    properties = {'!SIM.spectral.wave_min': .549, '!SIM.spectral.wave_max': .551, '!OBS.mjdobs': 60310.25, '!INST.include_diffuse_emissivity': False}
    commands = sim.UserCommands(use_instrument='ZShooter_v2', properties=properties)
    ignore = ['atmo_transmission', 'continuum_emission', 'airglow_and_interline_continuum', 'adc_residuals', 'shot_noise', 'readout_noise_selector']
    for document in commands.yaml_dicts:
        if 'effects' in document:
            document['effects'] = [effect for effect in document['effects'] if effect['name'] not in ignore]
    train = sim.OpticalTrain(commands)
    wave = np.linspace(5400., 5600., 4001)*u.AA
    spectrum = SourceSpectrum(Empirical1D, points=wave, lookup_table=.01*np.exp(-.5*((wave.value-5500)/1.)**2)*PHOTLAM)
    source = sim.Source(x=[0.], y=[.05], ref=[0], weight=[1.], spectra=[spectrum])
    for slit, seconds, redshift, weight in ((.7, 3600., 0., 1.), (.33, 300., .0003, 2.), (.7, 3600., 0., 1.)):
        source.spectra[0].z = redshift
        source.fields[0].field['weight'][0] = weight
        settings = {'!INST.vis_curr_slit': slit, '!INST.nir_curr_slit': slit, '!OBS.dit': seconds, '!OBS.ndit': 1}
        train.cmds.update(properties=settings.copy())
        train.update(migrate_cmds=True)
        fresh_commands = deepcopy(commands)
        fresh_commands.update(properties=settings)
        fresh = sim.OpticalTrain(fresh_commands)
        train.observe(source, update=False)
        fresh.observe(source, update=False)
        for actual, expected in zip(train.image_planes, fresh.image_planes):
            np.testing.assert_allclose(actual.data, expected.data, rtol=1e-6, atol=1e-12, err_msg=f'slit={slit}, DIT={seconds}')
        for arm in range(6):
            aperture = train['slitwheel_selector'].get_effect(arm).current_slit.table
            assert np.ptp(aperture['y']) == pytest.approx(slit)
        readouts = train.readout(dit=seconds, ndit=1, dit_green=seconds/2, ndit_green=2)
        for index, readout in enumerate(readouts):
            header = readout[1].header
            assert header['DIT'] == (seconds/2 if index == 1 else seconds)
            assert header['NDIT'] == (2 if index == 1 else 1)
            assert header['EXPTIME'] == seconds


@pytest.fixture(scope='module')
def spec_commands():
    with pytest.MonkeyPatch.context() as patch:
        patch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(ROOT.parent))
        yield sim.UserCommands(use_instrument='ZShooter_v2', set_modes=['SPEC'])


@pytest.fixture(scope='module')
def traces(spec_commands):
    config = yaml.safe_load((ROOT / 'ZShooter_SPEC.yaml').read_text())
    effect = next(item for item in config['effects'] if item['name'] == 'trace_list_analytical')
    kwargs = dict(effect['kwargs'])
    kwargs['filename'] = str(ROOT / kwargs['filename'])
    simulation = next(item for item in yaml.safe_load_all((ROOT / 'default.yaml').read_text()) if item['alias'] == 'SIM')
    kwargs['trace_flux_jacobian'] = simulation['properties']['spectral']['trace_flux_jacobian']
    return EchelleSpectralTraceList(cmds=spec_commands, **kwargs).spectral_traces


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


@pytest.mark.parametrize('detector_id', range(6))
def test_dark_current_participates_in_poisson_statistics(spec_commands, detector_id):
    config = yaml.safe_load((ROOT / 'ZShooter_DET.yaml').read_text())
    names = ['exposure_integration_selector', 'dark_current_selector', 'shot_noise', 'readout_noise_selector']
    effects = [item for item in config['effects'] if item['name'] in names]
    assert len(effects) == 4
    cmds = deepcopy(spec_commands)
    cmds['!OBS.dit'], cmds['!OBS.ndit'], cmds['!SIM.random.seed'] = 3600., 1, 126
    effects = [make_effect(item, cmds) for item in effects]
    dark_effect, = [effect.get_effect(detector_id) for effect in effects if effect.meta['name'] == 'dark_current_selector']
    read_effect, = [effect.get_effect(detector_id) for effect in effects if effect.meta['name'] == 'readout_noise_selector']
    dark = dark_effect._get_dark_level(detector_id)
    ron = from_currsys(read_effect.meta['noise_std'], cmds)
    header = header_from_list_of_xy([-4.5, 4.5], [-4.5, 4.5], .015, 'D')
    detector = Detector(header, cmds=cmds, id=detector_id)
    for effect in effects:
        effect.apply_to(detector)
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
def test_native_trace_photon_normalization_is_independent_of_sky_grid(traces, spec_commands, channel):
    trace = next(trace for name, trace in traces.items() if name.startswith(channel + '_'))
    center = (trace.wave_min + trace.wave_max) / 2
    # Native detector geometry comes from the same IRDB prescription as the trace.
    config = yaml.safe_load((ROOT / 'ZShooter_DET.yaml').read_text())
    detector_effect = next(item for item in config['effects'] if item['class'] == 'DetectorList' and item['kwargs']['image_plane_id'] == trace.meta['image_plane_id'])
    detector_header = make_effect(detector_effect, spec_commands).image_plane_header
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


@pytest.mark.parametrize('channel,center', [('b', .365), ('g', .53), ('r', .75), ('yj', 1.2), ('h', 1.65), ('k', 2.2)])
def test_observed_point_flux_converges_with_spatial_grid(monkeypatch, channel, center):
    from synphot import SourceSpectrum, ConstFlux1D
    from synphot.units import PHOTLAM

    monkeypatch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(ROOT.parent))
    # Keep the prescribed optics, effective PSF and native detector pitch.
    # Sky and ADC are excluded to isolate spatial integration rather than
    # changing a quantized atmospheric displacement between grids.
    ignore = ['atmo_transmission', 'continuum_emission', 'airglow_and_interline_continuum', 'adc_residuals']
    source = sim.Source(x=[.037], y=[.13], ref=[0], weight=[1.], spectra=[SourceSpectrum(ConstFlux1D, amplitude=.01*PHOTLAM)])
    rates = []
    sky_steps = []
    for step in [.1570545, .07852725, .039263625]:
        properties = {'!INST.pixel_scale': step, '!SIM.spectral.wave_min': center*.999, '!SIM.spectral.wave_max': center*1.001, '!OBS.seeing': .7, '!OBS.mjdobs': 60310.25, '!INST.vis_curr_slit': .7, '!INST.nir_curr_slit': .7, '!INST.include_diffuse_emissivity': False}
        properties.update({f'!DET.pixel_scale_{arm}': step for arm in ('b', 'g', 'r', 'yj', 'h', 'k')})
        commands = sim.UserCommands(use_instrument='ZShooter_v2', properties=properties)
        for document in commands.yaml_dicts:
            if 'effects' in document:
                document['effects'] = [effect for effect in document['effects'] if effect['name'] not in ignore]
        train = sim.OpticalTrain(commands)
        sky_steps.append(abs(train.fov_manager.fovs[0].header['CDELT1'])*3600)
        assert all(manager[0]._native_header['CDELT1D'] == .015 for manager in train.detector_managers)
        train.observe(None)
        sky = sum(plane.data.sum() for plane in train.image_planes)
        train.observe(source)
        rates.append(sum(plane.data.sum() for plane in train.image_planes) - sky)
    print(f'{channel}: source e/s {rates}; relative to finest {np.array(rates)/rates[-1]}')
    assert sky_steps[0] > sky_steps[1] > sky_steps[2]
    assert np.all(np.array(rates) > 0)
    np.testing.assert_allclose(np.array(rates)/rates[-1], 1., rtol=.02)
