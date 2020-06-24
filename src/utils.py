import _MEWtools as mew
import yaml
import pathlib
import os

STRUCTURE_PATH = 'structures'
INTERIOR_PATH = 'interiors'


def import_structure(name, overrides={}):
    folder = pathlib.Path(__file__).parent.absolute()
    path = os.path.join(folder, STRUCTURE_PATH, f'{name}.yml')

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

    g = 6.67e-11

    mass = float(config['mass']) if config.get('mass') is not None else 0
    semi_major_distance = float(config['semiMajorDistance']) if config.get('semiMajorDistance') is not None else 0
    eccentricity = float(config['eccentricity']) if config.get('eccentricity') is not None else 0
    obliquity = float(config['obliquity']) if config.get('obliquity') is not None else 0
    obliquity_phase = float(config['obliquityPhase']) if config.get('obliquityPhase') is not None else 0
    obliquity_phase_rate = float(config['obliquityPhaseRate']) if config.get('obliquityPhaseRate') is not None else 0

    spin_rate = float(config['spinRate']) if config.get('spinRate') is not None else 0

    non_synchronous_rotation_rate = float(config['nonSynchronousRotationRate']) if config.get(
        'nonSynchronousRotationRate') is not None else 0
    libration_amplitude = float(config['librationAmplitude']) if config.get('librationAmplitude') is not None else 0
    libration_phase = float(config['librationPhase']) if config.get('librationPhase') is not None else 0
    libration_frequency = float(config['librationFrequency']) if config.get('librationFrequency') is not None else 0

    sat = mew.satellite(
        rin=radii,
        pin=density,
        uin=rigidity,
        nin=viscosity,
        G=g,
        M=mass,
        a=semi_major_distance,
        e=eccentricity,
        o=obliquity,
        op=obliquity_phase,
        on=obliquity_phase_rate,
        Ω=spin_rate,
        NSRn=non_synchronous_rotation_rate,
        libA=libration_amplitude,
        libp=libration_phase,
        libn=libration_frequency
    )
    return sat


def clear_extension(file):
    return file.replace('.yml', '').replace('.yaml', '')


def list_config_files(config_path):
    folder = pathlib.Path(__file__).parent.absolute()
    path = os.path.join(folder, config_path)

    files = os.listdir(path)
    configs = list(map(clear_extension, files))

    for config in configs:
        print(config)

    return configs


def list_structures():
    return list_config_files(STRUCTURE_PATH)


def import_interior(name):
    folder = pathlib.Path(__file__).parent.absolute()
    path = os.path.join(folder, INTERIOR_PATH, f'{name}.yml')

    with open(path, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except:
            config = {}

    return Interior(config)


def list_interiors():
    return list_config_files(INTERIOR_PATH)


class Interior:

    def __init__(self, config):
        k_year = 1. / 3.156E10

        self.viscosity = config['viscosity']
        self.rigidity = config['rigidity']
        self.surface_gravity = config['surface_gravity']
        self.radius = config['radius']
        self.eccentricity = config['eccentricity']
        self.he = config['love_h_elastic']
        self.le = config['love_e_elastic']

        def calc_strength(value):
            return k_year * value

        self.modal_strengths = list(map(calc_strength, config['modal_strengths']))
        self.hv = config['love_h_viscoelastic']
        self.lv = config['love_l_viscoelastic']

        self.he_NSR = config.get('love_h_elastic_nsr')
        self.le_NSR = config.get('love_l_elastic_nsr')

        def invert(value):
            return 1.0 / value

        self.modal_strengths_nsr = list(map(invert, config['modal_strengths_nsr'])) \
            if config.get('modal_strengths_nsr') is not None \
            else []

        self.hv_NSR = config['love_h_viscoelastic_nsr'] \
            if config.get('love_h_viscoelastic_nsr') is not None \
            else []

        self.lv_NSR = config['love_l_viscoelastic_nsr'] \
            if config.get('love_l_viscoelastic_nsr') is not None \
            else []


def round_heading(value, base=5):
    return base * round(value / base)
