"""Prepare ScopeSim tables from the calibrated ESO X-shooter optical model.

All inputs are pinned in calibration/. The attributed ESO mathematical kernel
is compiled only for preparation. Simulation needs the generated FITS tables,
not an ESO installation or a compiler. Run this script from any directory.
"""
from pathlib import Path
import argparse
import tempfile

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from scopesim.effects.spectral_trace_list_utils import SpectralTrace
from calibration.physical_model import PhysicalModel, build_kernel

ROOT = Path(__file__).parent
# Full standard readout science geometry (overscan omitted for CCDs).
ARM = {'uvb': (2048, 3000, .015), 'vis': (2048, 4000, .015), 'nir': (2048, 1100, .018)}
# Gaussian optical blur fitted to Table 13 after slit and pixel integration.
OPTICAL_R = {'uvb': 17000, 'vis': 43000, 'nir': 46000}
# User Manual v10, Table 14, NIR orders 11--26; wavelengths in nm.
NIR_BLAZE = [2377.28, 2179.17, 2011.54, 1867.86, 1743.33, 1634.38, 1538.23, 1452.78, 1376.31, 1307.5, 1245.24, 1188.64, 1136.96, 1089.58, 1046.0, 1005.77]


def focal_plane(arm, xp, yp):
    """ESO internal half-pixel-centered coordinates -> raw-oriented focal mm.

    CCDs: undo pipeline UVB rotation, with prescan/overscan already removed.
    NIR: undo CPL's positive quarter-turn and the crop of raw columns 4:2044
    and rows 20:1076 (xsh_data_pre.c). Standard raw readout is 2048 x 1100.
    """
    nx, ny, pixel = ARM[arm]
    if arm == 'uvb':
        x, y = nx - xp, ny - yp
    elif arm == 'vis':
        x, y = xp, yp
    else:
        x, y = 2044 - yp, xp + 20
    return (x - nx / 2) * pixel, (y - ny / 2) * pixel


def prepare(build):
    library = build_kernel(build)
    optical_response = fits.open(ROOT / 'calibration' / 'optical_response.fits')
    catalog, traces, efficiencies, checks = [], [], [], []
    cache = ROOT / 'calibration' / 'eso_etc'
    for arm_id, arm in enumerate(ARM):
        nx, ny, pixel = ARM[arm]
        model = PhysicalModel(library, arm)
        orders = {'uvb': range(13, 25), 'vis': range(16, 31), 'nir': range(11, 27)}[arm]
        if arm == 'nir':
            optics = np.loadtxt(cache / 'instr_NIR.dat')
            correction = np.loadtxt(cache / 'NIR_fudge.dat')
        qe = np.loadtxt(cache / f'ccd{926 + arm_id}.dat')
        for index, order in enumerate(orders):
            blaze_nm = 2 * np.sin(np.deg2rad(-model.parameters['nug'])) / (order * model.parameters['sg']) * 1e6
            scan = np.linspace(blaze_nm * (1 - 2 / order), blaze_nm * (1 + 2 / order), 2049)
            covered = []
            for s in [-5.5, 0, 5.5]:
                xp, yp = model(scan, np.full_like(scan, s), order)
                x, y = focal_plane(arm, xp, yp)
                inside = (abs(x) <= nx * pixel / 2) & (abs(y) <= ny * pixel / 2)
                covered.extend(scan[inside])
            # One scan interval beyond the detector intersection keeps edge pixels
            # covered. The physical model, rather than an extrapolated inverse,
            # provides the coordinates outside the illuminated detector footprint.
            limits = [min(covered) - (scan[1] - scan[0]), max(covered) + (scan[1] - scan[0])]
            wave, slit = np.meshgrid(np.linspace(*limits, 193), np.linspace(-5.5, 5.5, 21))
            xp, yp = model(wave.ravel(), slit.ravel(), order)
            x, y = focal_plane(arm, xp, yp)
            table = Table([wave.ravel() / 1000 * u.um, slit.ravel() * u.arcsec, x * u.mm, y * u.mm], names=['wavelength', 's', 'x', 'y'])
            name = f'{arm.upper()}_{order}'
            hdu = fits.table_to_hdu(table)
            hdu.name = name
            hdu.header['DISPDIR'] = 'x' if arm == 'nir' else 'y'
            hdu.header['PIXSIZE'] = (pixel, 'Physical unbinned pixel size [mm]')
            wv, sv = np.meshgrid(np.linspace(*limits, 90)[1:-1], np.linspace(-5.4, 5.4, 12))
            xv, yv = model(wv.ravel(), sv.ravel(), order)
            xv, yv = focal_plane(arm, xv, yv)
            for degree in range(4, 13):
                hdu.header['FITDEG'] = degree
                trace = SpectralTrace(hdu, trace_flux_jacobian=True)
                dx = trace.xilam2x(sv.ravel(), wv.ravel() / 1000) - xv
                dy = trace.xilam2y(sv.ravel(), wv.ravel() / 1000) - yv
                residual = np.hypot(dx, dy) / pixel
                si, wi = trace.xy2xi(xv, yv), trace.xy2lam(xv, yv)
                inverse = np.hypot(trace.xilam2x(si, wi) - xv, trace.xilam2y(si, wi) - yv) / pixel
                if max(np.sqrt(np.mean(residual**2)), np.sqrt(np.mean(inverse**2))) < .03 and max(residual.max(), inverse.max()) < .1:
                    break
            else:
                raise ValueError(f'{name}: forward/inverse mapping fails the subpixel requirement')
            checks.append((name, degree, np.sqrt(np.mean(residual**2)), residual.max(), np.sqrt(np.mean(inverse**2)), inverse.max()))
            # The cross-slit wavelength shift comes from the same calibrated
            # optical model, including prism anamorphism. It is independent of
            # the numerical sky-grid spacing and the requested on-chip binning.
            xp, yp = model(wave.ravel(), slit.ravel(), order, eta_arcsec=.01)
            xm, ym = model(wave.ravel(), slit.ravel(), order, eta_arcsec=-.01)
            xp, yp = focal_plane(arm, xp, yp)
            xm, ym = focal_plane(arm, xm, ym)
            table['dlam_deta'] = (trace.xy2lam(xp, yp) - trace.xy2lam(xm, ym)) / .02 * u.um / u.arcsec
            header = hdu.header.copy()
            hdu = fits.table_to_hdu(table)
            hdu.header.extend(header, update=True)
            hdu.header['SLITCOL'] = 'dlam_deta'
            hdu.header['LSFRES'] = OPTICAL_R[arm]
            print(name, 'degree', degree, 'RMS [pix]', round(checks[-1][2], 5), flush=True)
            catalog.append((name, len(traces) + 2, arm_id, arm_id))
            traces.append(hdu)
            if arm == 'nir':
                # ETC envelope and scalar diffraction blaze are an explicit
                # approximation to the NIR order response.
                ew = np.linspace(*limits, 1025)
                blaze = np.sinc(order * (NIR_BLAZE[index] / ew - 1))**2
                throughput = np.interp(ew, optics[:, 0], optics[:, 1], left=0, right=0)
                throughput *= np.interp(ew, correction[:, 0], correction[:, 1]) * blaze
                eff = fits.table_to_hdu(Table([ew / 1000 * u.um, throughput], names=['wavelength', 'efficiency']))
                eff.name = name
            else:
                eff = optical_response[name].copy()
            efficiencies.append(eff)
        with (ROOT / 'throughput' / f'qe_{arm}.dat').open('w') as output:
            output.write('# wavelength_unit: um\n# Source: ESO ETC ccd' + str(926 + arm_id) + '.dat\nwavelength transmission\n')
            np.savetxt(output, np.column_stack([qe[:, 0] / 1000, qe[:, 1]]), fmt='%.9g')
    optical_response.close()
    cat = fits.table_to_hdu(Table(rows=catalog, names=['description', 'extension_id', 'aperture_id', 'image_plane_id']))
    primary = fits.PrimaryHDU()
    primary.header['ECAT'] = 1
    primary.header['EDATA'] = 2
    primary.header['INSTRUME'] = 'XSHOOTER'
    fits.HDUList([primary, cat] + traces).writeto(ROOT / 'traces' / 'traces.fits', overwrite=True, checksum=True)
    fits.HDUList([primary.copy(), cat.copy()] + efficiencies).writeto(ROOT / 'throughput' / 'orders.fits', overwrite=True, checksum=True)
    # Table 12 measures the warm stray light at four wavelengths, separately
    # from dark current and sky. Approximate its detector variation by a smooth
    # cross-dispersion profile, constant outside the measured positions. No
    # structure along dispersion is constrained by these four measurements.
    measurements = np.array([1.048, 1.238, 1.300, 1.682])
    positions = []
    nir_traces = [SpectralTrace(hdu) for hdu in traces if hdu.name.startswith('NIR')]
    for wave in measurements:
        trace = min(nir_traces, key=lambda tr: abs((tr.wave_min + tr.wave_max) / 2 - wave))
        positions.append(trace.xilam2y(0, wave) / ARM['nir'][2] + ARM['nir'][1] / 2 - .5)
    index = np.argsort(positions)
    positions = np.array(positions)[index]
    unfiltered = np.array([.06, .10, .13, .15])[index]
    filtered = np.array([.02, .03, .04, .04])[index]
    background = [fits.PrimaryHDU()]
    for slit_name in ['0.4', '0.6', '0.9', '1.2', '5.0', '0.6JH', '0.9JH']:
        aperture = Table.read(ROOT / 'slits' / f'slit_{slit_name}.dat', format='ascii')
        slit_width = np.ptp(aperture['y'])
        rates = filtered if slit_name.endswith('JH') else unfiltered
        row_rate = np.interp(np.arange(ARM['nir'][1]), positions, rates) * slit_width
        rate_map = np.broadcast_to(row_rate[:, None], (ARM['nir'][1], ARM['nir'][0])).astype('float32')
        hdu = fits.CompImageHDU(rate_map, name=slit_name, compression_type='GZIP_1', quantize_level=0)
        hdu.header['BUNIT'] = 'electron / s'
        background.append(hdu)
    fits.HDUList(background).writeto(ROOT / 'detector' / 'nir_background.fits', overwrite=True, checksum=True)
    Table(rows=checks, names=['order', 'degree', 'forward_rms_pix', 'forward_max_pix', 'inverse_rms_pix', 'inverse_max_pix']).write(ROOT / 'calibration' / 'trace_residuals.ecsv', overwrite=True)
    with (ROOT / 'detector' / 'nir_read_noise.dat').open('w') as output:
        output.write('# dit_unit: s\n# read_noise_unit: electron\n# Source: ESO ETC nir_dit_ron.dat\ndit read_noise\n')
        read_noise = np.loadtxt(cache / 'nir_dit_ron.dat')
        measured = (read_noise[:, 0] >= 1) & (read_noise[:, 0] <= 1800)
        np.savetxt(output, read_noise[measured], fmt='%.9g')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='xshooter-model-') as directory:
        prepare(Path(directory))
