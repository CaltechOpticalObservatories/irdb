import numpy as np
import astropy.units as u
from scipy.stats import norm
from scipy.signal import convolve
import yaml
import os
from synphot import (ThermalSpectralElement, Empirical1D)

def make_dichroic(
        domain: tuple[float | u.Quantity, float | u.Quantity],
        transmission: float,
        reflection: float,
        min: float | u.Quantity,
        max: float | u.Quantity,
        fold: float | u.Quantity,
        transition_width: float | u.Quantity = 30,
        temperature: float | u.Quantity = 0) -> tuple[ThermalSpectralElement, ThermalSpectralElement]:
    """
    Return a matched pait of dichroic efficiency curves.
    Args:
        domain:
        min:
        max:
        fold:
        transmission:
        reflection:
        transition_width:

    Returns: (ThermalSpectralElement, ThermalSpectralElement)

    """
    domain = u.Quantity(domain)
    full_x = np.linspace(np.floor(domain.min()),
                         np.ceil(domain.max()),
                         num=int(np.ceil(np.diff(domain)).value))
    full_yt = np.zeros(full_x.size)
    full_yr = np.zeros_like(full_yt)
    in_trans = (full_x >= fold)
    in_refl = (~in_trans)
    full_yt[in_trans & (full_x <= max)] = transmission
    full_yr[in_refl & (full_x >= min)] = reflection

    if not isinstance(transition_width, u.Quantity):
        transition_width = transition_width * domain.unit
    transition_sigma = transition_width / 2 / np.sqrt(2 * np.log(2))
    ker_f = norm(loc=0, scale=transition_sigma.value).pdf
    kx = np.arange(0, 10 * transition_width.value, (full_x[1] - full_x[0]).value)
    kx = np.concatenate([-kx[1:][::-1], kx])
    ker = ker_f(kx)

    # ker = ker_f(np.linspace(-10*transition_width, 10 * transition_width, num=n))
    smoothedt = convolve(full_yt, ker, mode='same')
    smoothedr = convolve(full_yr, ker, mode='same')

    return (ThermalSpectralElement(Empirical1D, points=full_x, lookup_table=smoothedr, fill_value=0, method='linear',
                                   bounds_error=False, temperature=temperature),
            ThermalSpectralElement(Empirical1D, points=full_x, lookup_table=smoothedt, fill_value=0, method='linear',
                                   bounds_error=False, temperature=temperature))


def make_zshooter_dichroics(channel_ranges:dict[str, tuple[float, float]], valid_domain: tuple[float, float],
                            transmission: float, reflection: float, range_pad: float):
    data = {'I': make_dichroic(valid_domain, transmission, reflection,
                                channel_ranges['U'].min() - range_pad, channel_ranges['K'].max() + range_pad,
                                0.5 * (channel_ranges['Y'].min() + channel_ranges['I'].max())),
             'U': make_dichroic(valid_domain, transmission, reflection,
                                channel_ranges['U'].min() - range_pad, channel_ranges['I'].max() + range_pad,
                                0.5 * (channel_ranges['B'].min() + channel_ranges['U'].max())),
             'G': make_dichroic(valid_domain, transmission, reflection,
                                channel_ranges['B'].min() - range_pad, channel_ranges['I'].max() + range_pad,
                                0.5 * (channel_ranges['R'].min() + channel_ranges['G'].max())),
             'B': make_dichroic(valid_domain, transmission, reflection,
                                channel_ranges['B'].min() - range_pad, channel_ranges['G'].max() + range_pad,
                                0.5 * (channel_ranges['G'].min() + channel_ranges['B'].max())),
             'R': make_dichroic(valid_domain, transmission, reflection,
                                channel_ranges['R'].min() - range_pad, channel_ranges['I'].max() + range_pad,
                                0.5 * (channel_ranges['I'].min() + channel_ranges['R'].max())),
             'H': make_dichroic(valid_domain, transmission, reflection,
                                channel_ranges['H'].min() - range_pad, channel_ranges['K'].max() + range_pad,
                                0.5 * (channel_ranges['K'].min() + channel_ranges['H'].max())),
             }
    if 'J' in channel_ranges:
        data.update({
            'Y': make_dichroic(valid_domain, transmission, reflection,
                               channel_ranges['Y'].min() - range_pad, channel_ranges['J'].max() + range_pad,
                               0.5 * (channel_ranges['J'].min() + channel_ranges['Y'].max())),
            'J': make_dichroic(valid_domain, transmission, reflection,
                               channel_ranges['Y'].min() - range_pad, channel_ranges['K'].max() + range_pad,
                               0.5 * (channel_ranges['H'].min() + channel_ranges['J'].max()))
                    })
    else:
        data['J']=make_dichroic(valid_domain, transmission, reflection,
                           channel_ranges['Y'].min() - range_pad, channel_ranges['K'].max() + range_pad,
                           0.5 * (channel_ranges['H'].min() + channel_ranges['Y'].max()))
    return data

DICHROIC_PAD = 30 * u.nm
DICHROIC_TRANS = 0.985
DICHROIC_REFL = 0.985
DICHROIC_DOMAIN = (250, 2550) * u.nm

# get path of this script, and use it to find the yaml config file in the same directory
channels_yaml_path = os.path.join(os.path.dirname(__file__), 'zshooter_channels.yaml')
zshooter_cfg = yaml.safe_load(open(channels_yaml_path, 'r'))['zshooter_channels']

channel_range = {}
for c, cfg in zshooter_cfg.items():
        crange = cfg['dispersion']['range'] * u.nm
        channel_range[c] = crange

dichroics = make_zshooter_dichroics(channel_range, DICHROIC_DOMAIN, DICHROIC_TRANS, DICHROIC_REFL, DICHROIC_PAD)