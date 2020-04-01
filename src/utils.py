import _MEWtools as MEW
import yaml
import pathlib
import os

def import_structure(name, overrides={}):
    folder = pathlib.Path(__file__).parent.absolute()
    path = os.path.join(folder, 'structures', f'{name}.yml')
    
    with open(path, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except:
            config = {}

    config.update(overrides)

    radii = list(map(float, config['radii']))
    density = list(map(float, config['density']))
    rigidity = list(map(float, config['rigidity']))
    viscosity = list(map(float, config['viscosity']))

    G = 6.67e-11

    mass = float(config['mass']) if config.get('mass') is not None else 0
    semiMajorDistance = float(config['semiMajorDistance']) if config.get('semiMajorDistance') is not None else 0
    eccentricity = float(config['eccentricity']) if config.get('eccentricity') is not None else 0
    obliquity = float(config['obliquity']) if config.get('obliquity') is not None else 0
    obliquityPhase = float(config['obliquityPhase']) if config.get('obliquityPhase') is not None else 0
    obliquityPhaseRate = float(config['obliquityPhaseRate']) if config.get('obliquityPhaseRate') is not None else 0

    spinRate = float(config['spinRate']) if config.get('spinRate') is not None else 0

    nonSynchronusRotationRate = float(config['nonSynchronusRotationRate']) if config.get('nonSynchronusRotationRate') is not None else 0
    librationAmplitude = float(config['librationAmplitude']) if config.get('librationAmplitude') is not None else 0
    librationPhase = float(config['librationPhase']) if config.get('librationPhase') is not None else 0
    librationFrequency = float(config['librationFrequency']) if config.get('librationFrequency') is not None else 0

    sat = MEW.satellite(
        rin = radii,
        pin = density,
        uin = rigidity,
        nin = viscosity,
        G = G,
        M = mass,
        a = semiMajorDistance,
        e = eccentricity,
        o = obliquity,
        op = obliquityPhase,
        on = obliquityPhaseRate,
        Ω = spinRate,
        NSRn = nonSynchronusRotationRate,
        libA = librationAmplitude,
        libp = librationPhase,
        libn = librationFrequency
    )
    return sat

def list_structures():
    folder = pathlib.Path(__file__).parent.absolute()
    path = os.path.join(folder, 'structures')

    files = os.listdir(path)
    clear_extension = lambda file: file.replace('.yml', '').replace('.yaml', '')
    structures = list(map(clear_extension, files))

    for structure in structures:
        print(structure)

    return structures
