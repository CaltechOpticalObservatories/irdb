"""Physical checks against the configured, full-frame X-shooter model."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

import scopesim as sim
from scopesim.detector import Detector
from scopesim.optics.image_plane_utils import header_from_list_of_xy

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def train(monkeypatch):
    monkeypatch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(ROOT.parent))
    cmds = sim.UserCommands(use_instrument='XSHOOTER', ignore_effects=['skycalc_atmosphere'])
    return sim.OpticalTrain(cmds)


def test_native_geometry_and_orders(train):
    assert [plane.data.shape for plane in train.image_planes] == [(3000, 2048), (4000, 2048), (1100, 2048)]
    traces = train['spectral_traces'].spectral_traces
    assert len(traces) == 43
    for arm, count, pitch in [('UVB', 12, .015), ('VIS', 15, .015), ('NIR', 16, .018)]:
        selected = [trace for name, trace in traces.items() if name.startswith(arm)]
        assert len(selected) == count
        assert all(trace.meta['pixel_size'] == pitch for trace in selected)
    assert train.cmds['!TEL.name'] == 'VLT UT3'


@pytest.mark.parametrize('bx,by', [(1, 1), (1, 2), (2, 2)])
def test_full_readout_charge_headers_and_repeat(train, bx, by):
    for arm in ['uvb', 'vis']:
        train.cmds[f'!DET.binx_{arm}'] = bx
        train.cmds[f'!DET.biny_{arm}'] = by
    for arm in ['uvb', 'vis', 'nir']:
        train.cmds[f'!OBS.dit_{arm}'] = 30.
        train.cmds[f'!OBS.ndit_{arm}'] = 3
    for name in ['shot_noise', 'readout_noise']:
        train[name].include = False
    for plane in train.image_planes:
        plane.data[:] = 2.
    for repeat in range(2):
        outputs = train.readout(reset=False)
        for index, arm in enumerate(['uvb', 'vis', 'nir']):
            hdu = outputs[index][1]
            header = hdu.header
            binx, biny = (bx, by) if arm != 'nir' else (1, 1)
            ny, nx = train.image_planes[index].data.shape
            gain = train.cmds[f'!DET.gain_{arm}']
            dark = train.cmds[f'!DET.dark_{arm}']
            integrations = 1 if arm == 'nir' else 3
            expected = (2 + dark)*30*integrations*binx*biny/gain
            assert hdu.data.shape == (ny//biny, nx//binx)
            np.testing.assert_allclose(hdu.data, expected, rtol=2e-6)
            assert outputs[index][0].header['ARM'] == arm.upper()
            assert header['DET_ID'] == index
            assert header['GAIN'] == gain
            assert header['BUNIT'] == 'adu'
            assert header['EXPTIME'] == 90
            assert header['COMBINE'] == ('average' if arm == 'nir' else 'sum')
            native = WCS(train.image_planes[index].header, key='D')
            actual = WCS(header, key='D').all_pix2world(13, 17, 0)
            np.testing.assert_allclose(actual, native.all_pix2world(13.5*binx-.5, 17.5*biny-.5, 0), atol=1e-10)


@pytest.mark.parametrize('arm,mode,gain,ron', [('uvb', 'slow', .63, 2.4), ('uvb', 'fast', 1.83, 4.6), ('vis', 'slow', .64, 3.4), ('vis', 'fast', 1.64, 6.6)])
def test_detector_modes(train, arm, mode, gain, ron):
    train.cmds[f'!DET.readout_mode_{arm}'] = mode
    train[f'{arm}_readout_parameters'].apply_to(None)
    assert train.cmds[f'!DET.gain_{arm}'] == gain
    assert train.cmds[f'!DET.ron_{arm}'] == ron


def test_all_slits_and_jh_background(train):
    for aid, arm in enumerate(['uvb', 'vis', 'nir']):
        wheel = train['slits'].wheel_effects[aid]
        for slit in wheel.meta['slit_names']:
            train.cmds[f'!INST.{arm}_curr_slit'] = slit
            aperture = wheel.current_slit.table
            assert np.ptp(aperture['x']) == 11
            assert np.ptp(aperture['y']) == pytest.approx(float(slit.removesuffix('JH')))
            if arm == 'nir':
                plane = train.image_planes[2]
                plane.data[:] = 0
                train['nir_stray_light'].apply_to(plane)
                width = np.ptp(aperture['y'])
                limits = (.02, .04) if slit.endswith('JH') else (.06, .15)
                assert plane.data.min()/width == pytest.approx(limits[0])
                assert plane.data.max()/width == pytest.approx(limits[1])
                curve = train['nir_order_blocking'].wheel_effects[2].current_filter.surface.transmission
                if slit.endswith('JH'):
                    assert curve(22000).value < .01
                    assert curve(15000).value == pytest.approx(.945)
                else:
                    assert curve(22000).value == 1


@pytest.mark.parametrize('aid,ndit', [(0, 1), (0, 4), (1, 4), (2, 1), (2, 4)])
def test_configured_noise_sequence(train, aid, ndit):
    arm = ['uvb', 'vis', 'nir'][aid]
    dit, rate = 100., .3
    train.cmds[f'!OBS.dit_{arm}'] = dit
    train.cmds[f'!OBS.ndit_{arm}'] = ndit
    train.cmds['!SIM.random.seed'] = 312
    if aid != 2:
        train[f'{arm}_readout_parameters'].apply_to(None)
        train.cmds[f'!DET.binx_{arm}'] = 2
        train.cmds[f'!DET.biny_{arm}'] = 2
    # Use a smaller detector only for statistical sample size; all effects and
    # arm parameters come from the actual full-frame train.
    header = header_from_list_of_xy([-6., 6.], [-4.5, 4.5], .015, 'D')
    det = Detector(header, cmds=train.cmds, id=aid)
    det.data[:] = rate
    for effect in train.optics_manager.detector_effects:
        effect.apply_to(det)
    gain, dark = train.cmds[f'!DET.gain_{arm}'], train.cmds[f'!DET.dark_{arm}']
    ron = det.header['RDNOISE'] if aid == 2 else train.cmds[f'!DET.ron_{arm}']
    area = 1 if aid == 2 else 4
    divisor = ndit if aid == 2 else 1
    mean = (rate + dark)*dit*ndit*area/divisor/gain
    variance = ((rate + dark)*dit*area + ron**2)*ndit/(divisor*gain)**2
    assert det.data.mean() == pytest.approx(mean, rel=.003)
    assert det.data.var() == pytest.approx(variance, rel=.025)


def test_unsupported_bin_fails(train):
    train.cmds['!DET.binx_uvb'], train.cmds['!DET.biny_uvb'] = 2, 1
    with pytest.raises(ValueError, match='Unsupported'):
        train.readout(reset=False)


def test_trace_fit_bounds_and_finite_physical_calibrations():
    residuals = Table.read(ROOT / 'calibration' / 'trace_residuals.ecsv')
    assert len(residuals) == 43
    for column in ['forward_rms_pix', 'inverse_rms_pix']:
        assert np.all(residuals[column] < .03)
    for column in ['forward_max_pix', 'inverse_max_pix']:
        assert np.all(residuals[column] < .1)
    with fits.open(ROOT / 'throughput' / 'orders.fits') as response:
        for hdu in response[2:]:
            assert np.all(np.isfinite(hdu.data['efficiency']))
            assert np.all(hdu.data['efficiency'] >= 0)
            assert np.all(hdu.data['efficiency'] <= 1)


def test_slit_resolving_power_against_manual():
    from XSHOOTER.calibration.fit_lsf import resolution_table
    resolution = resolution_table()
    assert np.all(abs(resolution['fractional_error']) < .08)
    for arm, width in [('UVB', 1.), ('VIS', .9), ('NIR', .9)]:
        selected = (resolution['arm'] == arm) & (resolution['slit_arcsec'] == width)
        assert abs(resolution['fractional_error'][selected][0]) < .03


@pytest.mark.parametrize('aid,name', [(0, 'UVB_18'), (1, 'VIS_23'), (2, 'NIR_16')])
def test_mapped_line_photons_centroid_and_variance(train, aid, name):
    trace = train['spectral_traces'].spectral_traces[name]
    center = (trace.wave_min + trace.wave_max)/2
    for step in [.1, .05]:
        wave = np.linspace(center*.999, center*1.001, 2001)
        xi = np.arange(-2, 2+step/2, step)
        eta = np.arange(-.5+step/2, .5, step)
        line = np.exp(-.5*((wave-center)/(center/100000))**2)
        illumination = np.exp(-.5*(eta/.2)**2)
        cube = line[:, None, None]*illumination[None, :, None]*np.exp(-.5*(xi/.35)**2)[None, None, :]
        header = fits.Header({'CTYPE1': 'LINEAR', 'CTYPE2': 'LINEAR', 'CTYPE3': 'WAVE', 'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec', 'CUNIT3': 'um', 'CDELT1': step, 'CDELT2': step, 'CDELT3': wave[1]-wave[0], 'CRPIX1': 1., 'CRPIX2': 1., 'CRPIX3': 1., 'CRVAL1': xi[0], 'CRVAL2': eta[0], 'CRVAL3': wave[0]})
        metadata = {'wave_min': wave[0]*u.um, 'wave_max': wave[-1]*u.um, 'xi_min': xi[0]*u.arcsec, 'xi_max': xi[-1]*u.arcsec, 'eta_min': -.5*u.arcsec, 'eta_max': .5*u.arcsec}
        fov = SimpleNamespace(cube=fits.ImageHDU(cube, header), header=header, detector_header=train.image_planes[aid].header, meta=metadata, trace_id=name)
        mapped = trace.map_spectra_to_focal_plane(fov)
        # Integrate source radiance over wavelength and both sky coordinates.
        photons = cube.sum()*step**2*(wave[1]-wave[0])
        assert mapped.data.sum() == pytest.approx(photons, rel=.005)
        yy, xx = np.indices(mapped.data.shape)
        xmm, ymm = WCS(mapped.header, key='D').all_pix2world(xx, yy, 0)
        detector_wave = trace.xy2lam(xmm, ymm)
        mean = (mapped.data*detector_wave).sum()/mapped.data.sum()
        variance = (mapped.data*(detector_wave-mean)**2).sum()/mapped.data.sum()
        eta_variance = (eta**2*illumination).sum()/illumination.sum()
        optical_sigma = center/trace.meta['optical_resolving_power']/np.sqrt(8*np.log(2))
        expected_variance = (center/100000)**2 + optical_sigma**2 + eta_variance*trace.slit_dispersion(center)**2
        assert abs(mean-center) < abs(trace.dlam_per_pix(center))*.02
        assert variance == pytest.approx(expected_variance, rel=.02)


@pytest.mark.parametrize('aid,center', [(0, .438), (1, .711), (2, 1.634)])
@pytest.mark.parametrize('point_source', [False, True])
def test_end_to_end_photon_rate(monkeypatch, aid, center, point_source):
    from synphot import SourceSpectrum
    from synphot.models import ConstFlux1D
    from synphot.units import PHOTLAM
    from scopesim.source.source_templates import uniform_source
    from scipy.integrate import quad
    from scipy.special import gamma

    monkeypatch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(ROOT.parent))
    lower, upper = center*.9995, center*1.0005
    ignore = ['skycalc_atmosphere', 'atmospheric_dispersion', 'nir_stray_light']
    if not point_source:
        ignore.append('seeing_psf')
    cmds = sim.UserCommands(use_instrument='XSHOOTER', ignore_effects=ignore, properties={'!SIM.spectral.wave_min': lower, '!SIM.spectral.wave_max': upper})
    train = sim.OpticalTrain(cmds)
    spectrum = SourceSpectrum(ConstFlux1D, amplitude=.01*PHOTLAM)
    source = sim.Source(x=[0], y=[0], ref=[0], weight=[1], spectra=[spectrum]) if point_source else uniform_source(spectrum, extent=20)
    train.observe(source)
    wave = np.linspace(lower, upper, 2001)*u.um
    throughput = np.zeros(len(wave))
    for name, trace in train['spectral_traces'].spectral_traces.items():
        if trace.meta['image_plane_id'] != aid:
            continue
        covered = (wave.value >= trace.wave_min) & (wave.value <= trace.wave_max)
        throughput[covered] += train['order_throughput'].efficiencies[name].throughput(wave[covered]).value
    qe = train['detector_qe'].wheel_effects[aid].throughput(wave).value
    telescope = train.optics_manager.get_all(sim.effects.SurfaceList)[0].throughput(wave).value
    aperture = train['slits'].wheel_effects[aid].current_slit.table
    if point_source:
        psf = train['seeing_psf']
        beta = psf.alpha
        scale = psf.fwhm(center*u.um).to_value(u.arcsec)/(2*np.sqrt(2**(1/beta)-1))
        norm = gamma(beta-.5)/(scale*np.sqrt(np.pi)*gamma(beta-1))
        width = np.ptp(aperture['y'])
        # The marginal Moffat profile integrates over the along-slit coordinate;
        # the 11 arcsec slit loses negligible flux at its ends for this seeing.
        slit_weight = quad(lambda eta: norm*(1+(eta/scale)**2)**(.5-beta), -width/2, width/2)[0]
    else:
        slit_weight = np.ptp(aperture['x'])*np.ptp(aperture['y'])
    density = spectrum(wave).to(u.photon/u.s/u.m**2/u.um)
    rate = np.trapezoid((density*train.optics_manager.area).value*throughput*qe*telescope, wave.value)*slit_weight
    # Finite sky and detector cells straddle the aperture/band edges. The source
    # is bright enough that telescope thermal emission is negligible here.
    assert train.image_planes[aid].data.sum() == pytest.approx(rate, rel=.01 if point_source else .025)
