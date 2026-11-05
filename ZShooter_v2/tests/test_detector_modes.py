from pathlib import Path

import numpy as np
import pytest
import scopesim as sim


@pytest.mark.parametrize('mode,detector_type,pitch,size', [('SPEC', 'LmAPD', .015, 2048), ('SPEC_ALTDETECT', 'H2RG', .018, 2048), ('SPEC_ALTDETECT4K', 'H4RG', .010, 4096)])
def test_k_detector_mode_reaches_geometry_and_readout_metadata(monkeypatch, mode, detector_type, pitch, size):
    root = Path(__file__).resolve().parents[2]
    monkeypatch.setitem(sim.rc.__config__, '!SIM.file.local_packages_path', str(root))
    commands = sim.UserCommands(use_instrument='ZShooter_v2', set_modes=[mode], properties={'!OBS.mjdobs': 60310.25})
    ignore = ['atmo_transmission', 'continuum_emission', 'airglow_and_interline_continuum', 'adc_residuals']
    for document in commands.yaml_dicts:
        if 'effects' in document:
            document['effects'] = [effect for effect in document['effects'] if effect['name'] not in ignore]
    train = sim.OpticalTrain(commands)
    assert train.cmds['!DET.detector_mode'] == mode
    detector, = train.detector_managers[5]
    header = detector._native_header
    assert header['DETMODE'] == mode
    assert header['DETTYPE'] == detector_type
    assert (header['NAXIS1'], header['NAXIS2']) == (size, size)
    np.testing.assert_allclose(abs(header['CDELT1D']), pitch)
    np.testing.assert_allclose(abs(header['CDELT2D']), pitch)
    traces = train['trace_list_analytical'].spectral_traces
    assert len([trace for trace in traces.values() if trace.meta['image_plane_id'] == 5]) == 8
