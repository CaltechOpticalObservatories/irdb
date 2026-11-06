"""Calibrate optical order response using the 2024-01-01 EG21 nod pair.

Run after prepare_calibrations.py, with the raw/reference files listed in
README.md in --data-dir. The first nod calibrates throughput; the second nod
checks repeatability. NIR measurements are retained for validation but are not
used to adjust the ETC response because of telluric and detector residuals.
"""
from pathlib import Path
import argparse
import numpy as np
from astropy import units as u, constants as c
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d, PchipInterpolator
from scopesim.effects.spectral_trace_list_utils import SpectralTrace

ROOT = Path(__file__).parent.parent


def order_rate(trace, difference, slit_position, pixel, gain, exposure, arm):
    """Sum a 3.6 arcsec spatial aperture in native detector rows or columns."""
    ny, nx = difference.shape
    wavelengths = np.linspace(trace.wave_min, trace.wave_max, 5000)
    slit = np.full_like(wavelengths, slit_position)
    x = trace.xilam2x(slit, wavelengths)/pixel + nx/2 - .5
    y = trace.xilam2y(slit, wavelengths)/pixel + ny/2 - .5
    dispersion, cross = (x, y) if arm == 'nir' else (y, x)
    order = np.argsort(dispersion)
    dp = np.arange(nx if arm == 'nir' else ny)
    covered = (dp > dispersion.min()+1) & (dp < dispersion.max()-1)
    dp = dp[covered]
    wave = interp1d(dispersion[order], wavelengths[order])(dp)
    cp = np.interp(dp, dispersion[order], cross[order])
    cross_pixels = np.rint(cp).astype(int)[:, None] + np.arange(-30, 31)
    on_detector = (cross_pixels >= 0) & (cross_pixels < (ny if arm == 'nir' else nx))
    cross_pixels = np.clip(cross_pixels, 0, (ny if arm == 'nir' else nx)-1)
    xx, yy = (dp[:, None]+np.zeros_like(cross_pixels), cross_pixels) if arm == 'nir' else (cross_pixels, dp[:, None]+np.zeros_like(cross_pixels))
    xmm, ymm = (xx+.5-nx/2)*pixel, (yy+.5-ny/2)*pixel
    s = trace.xy2xi(xmm, ymm)
    gx, gy = trace.xy2xi.gradient()
    ds = abs((gy if arm == 'nir' else gx)(xmm, ymm))*pixel
    weight = np.clip((1.8-abs(s-slit_position))/ds+.5, 0, 1)*on_detector
    rate = np.sum(difference[yy, xx]*weight, axis=1)*gain/exposure
    return wave, rate, abs(np.gradient(wave))


def measure(data_dir):
    traces = fits.open(ROOT/'traces/traces.fits')
    extinction = fits.getdata(data_dir/'atmos-ext-uvb.fits', 1)
    with fits.open(data_dir/'sky-default-raw.fits') as sky_file:
        sky = sky_file[1].data
        sky_wave = (sky['lam']*u.Unit(sky_file[1].columns['lam'].unit)).to_value(u.um)
    mirrors = Table.read(ROOT.parent/'VLT/TC_mirror_aluminium.dat', format='ascii')
    area = np.pi/4*(8.2**2-1.116**2)*u.m**2
    corrections, measurements, checks = [fits.PrimaryHDU()], [fits.PrimaryHDU()], []
    for arm_id, arm in enumerate(['uvb', 'vis', 'nir']):
        first = fits.open(data_dir/f'standard-{arm}.fits')
        second = fits.open(data_dir/f'standard-b-{arm}.fits')
        data_a, data_b = first[0].data, second[0].data
        if arm != 'nir':
            start = first[0].header['ESO DET OUT1 PRSCX']
            data_a, data_b = data_a[:, start:start+2048], data_b[:, start:start+2048]
        difference = data_a.astype('float64')-data_b
        reference = fits.getdata(data_dir/f'response-{arm}.fits', 1)
        qe = Table.read(ROOT/f'throughput/qe_{arm}.dat', format='ascii')
        for row in traces[1].data:
            if row['aperture_id'] != arm_id:
                continue
            trace = SpectralTrace(traces[row['extension_id']])
            measured = []
            for slit_position, sign, header in [(-2.5, 1, first[0].header), (2.5, -1, second[0].header)]:
                wave, rate, dwave = order_rate(trace, sign*difference, slit_position, trace.meta['pixel_size'], [.63, .64, 2.29][arm_id], header['EXPTIME'], arm)
                use = (wave*1000 >= reference['LAMBDA'][0]) & (wave*1000 <= reference['LAMBDA'][-1])
                wave, rate, dwave = wave[use], rate[use], dwave[use]
                flux = np.interp(wave*1000, reference['LAMBDA'], reference['REF'])*u.erg/u.cm**2/u.s/u.AA
                incident = (flux*area*(dwave*u.um)/(c.h*c.c/(wave*u.um))).to_value(1/u.s)
                airmass = header['ESO TEL AIRM START']
                if arm == 'nir':
                    transmission = np.interp(wave, sky_wave, sky['trans'])**(airmass/1.3)
                else:
                    transmission = 10**(-.4*airmass*np.interp(wave*1000, extinction['LAMBDA'], extinction['EXTINCTION']))
                telescope = np.interp(wave, mirrors['wavelength'], mirrors['reflection'])**2
                quantum_efficiency = np.interp(wave, qe['wavelength'], qe['transmission'])
                throughput = rate/(incident*transmission*telescope*quantum_efficiency)
                order = np.argsort(wave)
                measured.append((wave[order], throughput[order]))
            wave, eta_a = measured[0]
            wave_b, eta_b = measured[1]
            overlap = (wave >= wave_b[0]) & (wave <= wave_b[-1])
            wave, eta_a = wave[overlap], eta_a[overlap]
            eta_b = np.interp(wave, wave_b, eta_b)
            # Exclude molecular absorption and use robust 32-pixel continuum bins
            # to suppress cosmic rays and fixed-pattern structure in unflattened data.
            molecular = np.interp(wave, sky_wave, sky['trans_ma'])
            rows = []
            for start in range(0, len(wave)-31, 32):
                index = slice(start, start+32)
                if np.min(molecular[index]) < .97 or np.min(eta_a[index]) <= 0 or np.min(eta_b[index]) <= 0:
                    continue
                rows.append((np.mean(wave[index]), np.median(eta_a[index]), np.median(eta_b[index])))
            table = Table(rows=rows, names=['wavelength', 'throughput_a', 'throughput_b'])
            table['wavelength'].unit = u.um
            measurements.append(fits.table_to_hdu(table))
            measurements[-1].name = trace.trace_id
            if arm == 'nir':
                continue
            # Calibrate the smooth order response directly in physical throughput.
            # Outside the sampled interval, retain its endpoint value explicitly.
            # Molecular gaps are interpolated; they are not independent calibrations.
            profile = PchipInterpolator(table['wavelength'], table['throughput_a'])
            ew = np.linspace(trace.wave_min, trace.wave_max, 1025)
            efficiency = profile(np.clip(ew, table['wavelength'][0], table['wavelength'][-1]))
            output = fits.table_to_hdu(Table([ew*u.um, efficiency], names=['wavelength', 'efficiency']))
            output.name = trace.trace_id
            corrections.append(output)
            residual = table['throughput_b']/profile(table['wavelength'])-1
            checks.append((trace.trace_id, len(table), np.median(residual), np.sqrt(np.mean(residual**2))))
    fits.HDUList(corrections).writeto(ROOT/'calibration/optical_response.fits', overwrite=True, checksum=True)
    fits.HDUList(measurements).writeto(ROOT/'calibration/standard_measurements.fits', overwrite=True, checksum=True)
    Table(rows=checks, names=['order', 'continuum_bins', 'nod_b_fractional_bias', 'nod_b_fractional_rms']).write(ROOT/'calibration/standard_repeatability.ecsv', overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', required=True, type=Path)
    measure(parser.parse_args().data_dir)
