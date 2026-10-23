"""ESO's calibrated optical kernel, used only to prepare ScopeSim trace tables.

The library holds one active arm at a time, matching the preparation loop.
Simulation uses the prepared FITS files and requires neither this module nor a
C compiler. Internal model coordinates are converted separately by the caller.
"""
from pathlib import Path
import ctypes
import subprocess
import numpy as np
from astropy.table import Table

ROOT = Path(__file__).parent


def build_kernel(directory):
    directory.mkdir(parents=True, exist_ok=True)
    library = directory / 'eso_model.so'
    subprocess.run(['cc', '-shared', '-fPIC', '-O2', str(ROOT / 'eso_model.c'), '-o', str(library), '-lm'], check=True)
    return library


class PhysicalModel:
    def __init__(self, library, arm):
        self.library = ctypes.CDLL(str(library))
        self.library.set_parameter.argtypes = [ctypes.c_char_p, ctypes.c_double]
        self.library.reset_model.restype = None
        self.library.initialize_model.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        self.library.initialize_model.restype = None
        array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        self.library.evaluate.argtypes = [ctypes.c_int, array, array, ctypes.c_double, ctypes.c_int, array, array]
        self.library.evaluate.restype = None
        parameters = Table.read(ROOT / f'model_{arm}.fits', hdu=1)
        self.parameters = dict(zip(np.char.strip(parameters['Parameter_Name']), parameters['Best_Guess']))
        self.library.reset_model()
        angular = {'mues', 'nues', 'taues', 'cmup1', 'mug', 'nug', 'taug', 'mud', 'nud', 'taud', 'chiprot'}
        angular.update(f'{prefix}{i}' for prefix in ['mup', 'nup', 'taup'] for i in range(1, 7))
        for name, value in self.parameters.items():
            if name in angular:
                value = np.deg2rad(value)
            if self.library.set_parameter(name.encode(), value) != 1:
                raise KeyError(name)
        # xsh_data_instrument.c and xsh_model_arm_constants.h: internal model sizes.
        dimensions = {'uvb': (2048, 3000, 96., 0.), 'vis': (2048, 4000, 96., 0.), 'nir': (1020, 2040, 4., 8.)}
        self.library.initialize_model(*dimensions[arm])

    def __call__(self, wavelength_nm, slit_arcsec, order, eta_arcsec=0):
        wavelength_nm = np.ascontiguousarray(wavelength_nm, dtype=np.float64)
        slit_arcsec = np.ascontiguousarray(slit_arcsec, dtype=np.float64)
        x, y = np.empty_like(wavelength_nm), np.empty_like(wavelength_nm)
        self.library.evaluate(len(x), wavelength_nm, slit_arcsec, eta_arcsec, order, x, y)
        return x, y
