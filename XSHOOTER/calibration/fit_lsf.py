"""Fit optical Gaussian blur to the manual's slit-limited arc-line resolution.

Uniform slit illumination and the physical detector pixel are integrated before
fitting a Gaussian, as in an arc-line resolution measurement. Each arm uses the
mean slit dispersion and pixel sampling at its modeled order centers. This is
an arm-level approximation; it does not impose constant final resolving power.
"""
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit, least_squares
from scipy.special import erf
from scopesim.effects.spectral_trace_list_utils import SpectralTrace

ROOT = Path(__file__).resolve().parents[1]
# X-shooter User Manual v10, Table 13 (slit width [arcsec], resolving power).
REFERENCE = {
    'UVB': ([.5, .8, 1., 1.3, 1.6], [9700, 6700, 5400, 4100, 3200]),
    'VIS': ([.4, .7, .9, 1.2, 1.5], [18400, 11400, 8900, 6500, 5000]),
    'NIR': ([.4, .6, .9, 1.2], [11600, 8100, 5600, 4300]),
}


def gaussian(wave, amplitude, sigma):
    return amplitude*np.exp(-.5*(wave/sigma)**2)


def slit_line_resolving_power(width, r_slit, r_pixel, r_optical):
    # Fractional wavelength coordinate; both slit and detector pixel widths
    # are obtained from the configured optical model, not numerical sky pixels.
    sigma = 1/r_optical/np.sqrt(8*np.log(2))
    slit_width = width/r_slit
    wave = np.linspace(-3*slit_width-8*sigma, 3*slit_width+8*sigma, 1601)
    pixel = np.linspace(-.5/r_pixel, .5/r_pixel, 31)
    pos = wave[:, None] + pixel
    line = np.mean(erf((pos+slit_width/2)/(np.sqrt(2)*sigma)) - erf((pos-slit_width/2)/(np.sqrt(2)*sigma)), axis=1)
    parameters, _ = curve_fit(gaussian, wave, line, p0=[line.max(), slit_width/2.355])
    return 1/(abs(parameters[1])*np.sqrt(8*np.log(2)))


def resolution_table():
    rows = []
    with fits.open(ROOT / 'traces' / 'traces.fits') as hdul:
        for arm, (widths, reference) in REFERENCE.items():
            slit_sampling, pixel_sampling = [], []
            traces = [SpectralTrace(hdu) for hdu in hdul[2:] if hdu.name.startswith(arm)]
            for trace in traces:
                wave = (trace.wave_min + trace.wave_max)/2
                x, y = trace.xilam2x(0, wave), trace.xilam2y(0, wave)
                gx, gy = trace.xy2lam.gradient()
                slit_sampling.append(wave/abs(trace.slit_dispersion(wave)))
                pixel_sampling.append(wave/(np.hypot(gx(x, y), gy(x, y))*trace.meta['pixel_size']))
            r_slit, r_pixel = np.mean(slit_sampling), np.mean(pixel_sampling)
            configured_r = traces[0].meta['optical_resolving_power']

            def residual(log_r):
                modeled = [slit_line_resolving_power(width, r_slit, r_pixel, np.exp(log_r[0])) for width in widths]
                return np.array(modeled)/reference - 1

            fit = least_squares(residual, [np.log(configured_r)])
            fitted_r = np.exp(fit.x[0])
            print(f'{arm}: fitted optical R = {fitted_r:.0f}; configured R = {configured_r}')
            for width, published_r in zip(widths, reference):
                modeled_r = slit_line_resolving_power(width, r_slit, r_pixel, configured_r)
                rows.append((arm, width, published_r, modeled_r, modeled_r/published_r - 1))
    table = Table(rows=rows, names=['arm', 'slit_arcsec', 'reference_R', 'model_R', 'fractional_error'])
    table['slit_arcsec'].unit = 'arcsec'
    return table


if __name__ == '__main__':
    table = resolution_table()
    table.write(ROOT / 'calibration' / 'slit_resolution.ecsv', overwrite=True)
    print(table)
