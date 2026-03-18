import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import os


"""
Create a TRACE_ZSHOOTER.fits file from text files defining echelle orders
This only creates the trace definition file, no efficiency curves
"""
basepath = os.path.dirname(os.path.abspath(__file__))
f = fits.open(f'{basepath}/echelle_traces.fits')
hdr0 = f[0].header
hdr1 = f[1].header
tmphdr = f[2].header

ub_bottom = pd.read_csv(f'{basepath}/Echelle_Layout/Blue_Echelle_Bottom.txt', sep='\s+', header=None, skiprows=1)
ub_center = pd.read_csv(f'{basepath}/Echelle_Layout/Blue_Echelle_Center.txt', sep='\s+', header=None, skiprows=1)
ub_top = pd.read_csv(f'{basepath}/Echelle_Layout/Blue_Echelle_Top.txt', sep='\s+', header=None, skiprows=1)
vis_bottom = pd.read_csv(f'{basepath}/Echelle_Layout/Red_Echelle_Bottom.txt', sep='\s+', header=None, skiprows=1)
vis_center = pd.read_csv(f'{basepath}/Echelle_Layout/Red_Echelle_Center.txt', sep='\s+', header=None, skiprows=1)
vis_top = pd.read_csv(f'{basepath}/Echelle_Layout/Red_Echelle_Top.txt', sep='\s+', header=None, skiprows=1)
hk_bottom = pd.read_csv(f'{basepath}/Echelle_Layout/HK_Echelle_Bottom.txt', sep='\s+', header=None, skiprows=1)
hk_center = pd.read_csv(f'{basepath}/Echelle_Layout/HK_Echelle_Center.txt', sep='\s+', header=None, skiprows=1)
hk_top = pd.read_csv(f'{basepath}/Echelle_Layout/HK_Echelle_Top.txt', sep='\s+', header=None, skiprows=1)
yj_bottom = pd.read_csv(f'{basepath}/Echelle_Layout/YJ_Echelle_Bottom.txt', sep='\s+', header=None, skiprows=1)
yj_center = pd.read_csv(f'{basepath}/Echelle_Layout/YJ_Echelle_Center.txt', sep='\s+', header=None, skiprows=1)
yj_top = pd.read_csv(f'{basepath}/Echelle_Layout/YJ_Echelle_Top.txt', sep='\s+', header=None, skiprows=1)

hdulist = []
tracelist = []
extid = 2
# NIR HK
hk_orders = hk_bottom[0].unique()
for hko in hk_orders:
    wavelengths, xis, xs, ys = [], [], [], []
    for i,df in enumerate([hk_bottom, hk_center, hk_top]):
        subdf = df[df[0]==hko].reset_index(drop=True)
        subdf['xi'] = -5.0 + i*5.0
        subdf[1] = subdf[1]/1000.0  # convert to um
        subdf[2] = subdf[2] + 30.720  # shift to match detector coords
        wavelengths.extend(subdf[1].values)
        xis.extend(subdf['xi'].values)
        xs.extend(subdf[2].values)
        ys.extend(subdf[3].values)
    tab = Table({'wavelength': wavelengths*u.um,'s': xis*u.arcsec,'x': xs*u.mm,'y': ys*u.mm})
    hdr = tmphdr.copy()
    hdr['EXTNAME'] = f'hk_{int(hko)}'
    hdulist.append(fits.BinTableHDU(data=tab, header=hdr))
    tracelist.append([hdr['EXTNAME'], extid, 0, 0])
    extid += 1
# NIR YJ
yj_orders = yj_bottom[0].unique()
for yjo in yj_orders:
    wavelengths, xis, xs, ys = [], [], [], []
    for i,df in enumerate([yj_bottom, yj_center, yj_top]):
        subdf = df[df[0]==yjo].reset_index(drop=True)
        subdf['xi'] = -5.0 + i*5.0
        subdf[1] = subdf[1]/1000.0  # convert to um
        subdf[2] = subdf[2] - 30.720 # shift to match detector coords
        wavelengths.extend(subdf[1].values)
        xis.extend(subdf['xi'].values)
        xs.extend(subdf[2].values)
        ys.extend(subdf[3].values)
    tab = Table({'wavelength': wavelengths*u.um,'s': xis*u.arcsec,'x': xs*u.mm,'y': ys*u.mm})
    hdr = tmphdr.copy()
    hdr['EXTNAME'] = f'yj_{int(yjo)}'
    hdulist.append(fits.BinTableHDU(data=tab, header=hdr))
    tracelist.append([hdr['EXTNAME'], extid, 0, 0])
    extid += 1
# VIS
vis_orders = vis_bottom[0].unique()
for viso in vis_orders:
    wavelengths, xis, xs, ys = [], [], [], []
    for i,df in enumerate([vis_bottom, vis_center, vis_top]):
        subdf = df[df[0]==viso].reset_index(drop=True)
        subdf['xi'] = -5.0 + i*5.0
        subdf[1] = subdf[1]/1000.0  # convert
        wavelengths.extend(subdf[1].values)
        xis.extend(subdf['xi'].values)
        xs.extend(subdf[2].values)
        ys.extend(subdf[3].values)
    tab = Table({'wavelength': wavelengths*u.um,'s': xis*u.arcsec,'x': xs*u.mm,'y': ys*u.mm})
    hdr = tmphdr.copy()
    hdr['EXTNAME'] = f'gri_{int(viso)}'
    hdulist.append(fits.BinTableHDU(data=tab, header=hdr))
    tracelist.append([hdr['EXTNAME'], extid, 1, 1])
    extid += 1
# UB
ub_orders = ub_bottom[0].unique()
for ubo in ub_orders:
    wavelengths, xis, xs, ys = [], [], [], []
    for i,df in enumerate([ub_bottom, ub_center, ub_top]):
        subdf = df[df[0]==ubo].reset_index(drop=True)
        subdf['xi'] = -5.0 + i*5.0
        subdf[1] = subdf[1]/1000.0  # convert to um
        wavelengths.extend(subdf[1].values)
        xis.extend(subdf['xi'].values)
        xs.extend(subdf[2].values)
        ys.extend(subdf[3].values)
    tab = Table({'wavelength': wavelengths*u.um,'s': xis*u.arcsec,'x': xs*u.mm,'y': ys*u.mm})
    hdr = tmphdr.copy()
    hdr['EXTNAME'] = f'ub_{int(ubo)}'
    hdulist.append(fits.BinTableHDU(data=tab, header=hdr))
    tracelist.append([hdr['EXTNAME'], extid, 2, 2])
    extid += 1

hdul = fits.HDUList()
hdul.append(fits.PrimaryHDU(header=hdr0))
# Create trace list table
trace_tab = Table(rows=tracelist, names=['description', 'extension_id', 'aperture_id', 'image_plane_id'])
hdul.append(fits.BinTableHDU(data=trace_tab, header=hdr1))
for hdu in hdulist:
    hdul.append(hdu)

hdul.writeto(f'{basepath}/TRACE_ZSHOOTER.fits', overwrite=True)

"""
Create placeholder efficiency file TRACE_eff.fits
"""
hdul_eff = fits.HDUList()
hdul_eff.append(fits.PrimaryHDU(header=hdr0)) # same header as trace file
# create trace eff list table
trace_tab = np.array(tracelist)[:,:2].tolist()
trace_eff_tab = Table(rows=trace_tab, names=['description', 'extension_id'], dtype=['>S6', '>i4'])
trace_eff_tab['extension_id'] = trace_eff_tab['extension_id'].astype(int)
hdul_eff.append(fits.BinTableHDU(data=trace_eff_tab, header=hdr1))

for hdu in hdulist:
    tab = hdu.data
    wave = np.sort(np.unique(tab['wavelength']))
    # upside down quadratic efficiency curve, max efficiency 1.0 at center wavelength, 0 at edges
    eff = (wave - np.min(wave)) * (np.max(wave) - wave)
    eff = eff / np.max(eff)
    eff_tab = Table({'wavelength': wave*u.um, 'efficiency': eff})
    hdr = hdu.header.copy()
    hdul_eff.append(fits.BinTableHDU(data=eff_tab, header=hdr))
hdul_eff.writeto(f'{basepath}/TRACE_eff.fits', overwrite=True)
