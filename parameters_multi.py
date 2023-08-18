#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Memory Stimulation and Phase-Amplitude Coupling
# Copyright 2021 Nikolaos Vardlakis & Nicolas P. Rougier
# Released under the BSD 3-clauses license
# -----------------------------------------------------------------------------
import os
import json
import time
import subprocess
import numpy as np
from numpy import pi

# Constants
# noise_exc = 10. # uV
# noise_inh = 1. # uV
noise_exc = 0. # uV
noise_inh = 0. # uV
noise_EC = noise_DG = noise_CA3 = noise_CA1 = 0.

noise_EC_exc = noise_DG_exc = noise_CA3_exc = noise_CA1_exc = 0.0
noise_EC_inh = noise_DG_inh = noise_CA3_inh = noise_CA1_inh = 0.0
noise_EC_exc = 0. # 0.0031
noise_EC_inh = 0. # 0.0017
noise_DG_exc = 0. # 0.003
noise_DG_inh = 0. #
noise_CA3_exc = 0. # 0.003
noise_CA3_inh = 0. #
noise_CA1_exc = 0. # 0.0031
noise_CA1_inh = 0. #

# Default parameters
a = b = c = d = 0. # connections
a = 13. # 13
b = 0.14 # 0.14
c = 1.1 # 1.1
d = 0.2 # 0.2 // 0.43 new val
# a = 240. # 13
# b = 0.09 # 0.14
# c = 1.1 # 1.1
# d = 0.21 # 0.2 // 0.43 new val
I_in = 0.20 # 0.22 input gain Kuramoto
stim_amplitude = [9.] # nA
stim_onset = 1.8 # sec

# Default parameters
_data = {
    "seed_val"  : 42,       # Reproducibility

    # areas, tables 3.1-3.3, pages 45-48, Aussel
    "areas": {
        "EC"    : {
            "E" : {
                "N"     : int(10e3),
                "type"  : "PyCAN",
                "noise" : noise_exc     # Volts
            },
            "I" : {
                "N"     : int(1e3),
                "type"  : "Inh",
                "noise" : noise_inh
            }
        },
        "DG"    : {
            "E" : {
                "N"     : int(10e3),
                "type"  : "Py",
                "noise" : noise_exc
            },
            "I" : {
                "N"     : int(0.1e3),
                "type"  : "Inh",
                "noise" : noise_inh
            }
        },
        "CA3"   : {
            "E" : {
                "N"     : int(1e3),
                "type"  : "PyCAN",
                "noise" : noise_exc
            },
            "I" : {
                "N"     : int(0.1e3),
                "type"  : "Inh",
                "noise" : noise_inh
            }
        },
        "CA1"   : {
            "E" : {
                "N"     : int(10e3),
                "type"  : "PyCAN",
                "noise" : noise_exc
            },
            "I" : {
                "N"     : int(1e3),
                "type"  : "Inh",
                "noise" : noise_inh
            }
        }
    },

    # connectivity parameters
    "connectivity" : {
        "intra" : { # intra-area conn. probabilities per area |
            "EC"        : [[0., 0.37], [0.54, 0.]], # [[E-E, E-I], [I-E, I-I]]
            "DG"        : [[0., 0.06], [0.14, 0.]],
            "CA3"       : [[0.56, 0.75], [0.75, 0.]],
            "CA1"       : [[0., 0.28], [0.3, 0.7]]
        },
        "inter" : { # inter-area conn. probabilities
            "p_tri"     : 0.45,     # tri: [DG->CA3, CA3->CA1, CA1->EC] Aussel, pages 49,59
            "p_mono"    : 0.2       # mono: [EC->CA3, EC->CA1]
        },
        "inter_custom" : {
            "EC" : {
                "E" : [[0., 0.], [a, a], [b, b], [c, c]],
                "I" : [[0., 0.], [0., 0.], [0., 0.], [0., 0.]]
            },
            "DG" : {
                "E" : [[0., 0.], [0., 0.], [b, b], [0., 0.]],
                "I" : [[0., 0.], [0., 0.], [0., 0.], [0., 0.]]
            },
            "CA3" : {
                "E" : [[0., 0.], [0., 0.], [0., 0.], [c, c]],
                "I" : [[0., 0.], [0., 0.], [0., 0.], [0., 0.]]
            },
            "CA1" : {
                "E" : [[d, d], [0., 0.], [0., 0.], [0., 0.]],
                "I" : [[0., 0.], [0., 0.], [0., 0.], [0., 0.]]
            }
        }
    },

    # synapses
    # "synapses" : {
    #     "gmax_e" : 600.,    # pSiemens
    #     "gmax_i" : 60.
    # },

    # input parameters
    "fixed_input" : {
        "enabled"       : 0,                # [1=yes | 0=no]
        "low"           : 0.0,              # A0
        "high"          : 0.0,              # A1
        "frequency"     : 6.0,              # Hz
        "delay"         : 0.0,              # seconds
    },

    # Kuramoto oscillator parameters
    "Kuramoto" : {
        "N"             : 250,
        "f0"            : 6.,
        "sigma"         : 0.5,  # normal std
        "kN"            : 15,
        "gain_reset"    : 4.0,  # maybe 2?
        "gain_rhythm"   : np.around(I_in, 2), # nA
        "offset"        : -0*pi/2
    },

    # stimulation parameters
    "stimulation" : {
        "target"        : "CA1",            # target area [EC | DG | CA3 | CA1]
        "coordinates"   : (5.0, -8., 7.5),  # point electrode coordinates (x,y,z) [mm]
        "sigma"         : 0.33,             # conductivity of homogeneous conductive medium [S/m]
        "duration"      : 5.,               # [sec]
        "dt"            : .1e-3,            # [sec]
        "onset"         : stim_onset,       # [sec]
        "I"             : stim_amplitude,   # stimulation amplitude [nA]
        "pulse_width"   : [1.e-3],          # width (in time) of pulse ON phase [sec]
        "stim_freq"     : 5,                # stimulation frequency [Hz]
        "pulse_freq"    : 100,              # pulse frequency, determines ON duration [Hz]
        "nr_of_trains"  : 1,                # number of pulse trains
        "nr_of_pulses"  : 1,                # number of pulses per train
        "ipi"           : .1e-3             # inter-pulse interval [sec]
        },

    # simulation parameters
    "simulation" : {
        "duration"      : 10.0,             # second
        "dt"            : .1e-3,            # second
        "debugging"     : False
    },

    # git stuff
    "timestamp"         : None,
    "git_branch"        : None,
    "git_hash"          : None,
    "git_short_hash"    : None
}

def is_git_repo():
    """ Return whether current directory is a git directory """
    if subprocess.call(["git", "branch"],
            stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w')) != 0:
        return False
    return True

def get_git_revision_hash():
    """ Get current git hash """
    if is_git_repo():
        answer = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'])
        return answer.decode("utf8").strip("\n")
    return "None"

def get_git_revision_short_hash():
    """ Get current git short hash """
    if is_git_repo():
        answer = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'])
        return answer.decode("utf8").strip("\n")
    return "None"

def get_git_revision_branch():
    """ Get current git branch """
    if is_git_repo():
        answer = subprocess.check_output(
            ['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
        return answer.decode("utf8").strip("\n")
    return "None"

def default():
    """ Get default parameters """
    _data["timestamp"] = time.ctime()
    _data["git_branch"] = get_git_revision_branch()
    _data["git_hash"] = get_git_revision_hash()
    _data["git_short_hash"] = get_git_revision_short_hash()
    return _data

def save(filename, data=None):
    """ Save parameters into a json file """
    if data is None:
       data = { name : eval(name) for name in _data.keys()
                if name not in ["timestamp", "git_branch", "git_hash"] }
    data["timestamp"] = time.ctime()
    data["git_branch"] = get_git_revision_branch()
    data["git_hash"] = get_git_revision_hash()
    data["git_short_hash"] = get_git_revision_short_hash()
    with open(filename, "w") as outfile:
        json.dump(data, outfile, indent=4, sort_keys=False)

def load(filename):
    """ Load parameters from a json file """
    with open(filename) as infile:
        data = json.load(infile)
    return data

def dump(data):
    if not _data["timestamp"]:
        _data["timestamp"] = time.ctime()
    if not _data["git_branch"]:
        _data["git_branch"] = get_git_revision_branch()
    if not _data["git_hash"]:
        _data["git_hash"] = get_git_revision_hash()
        _data["git_short_hash"] = get_git_revision_short_hash()
    for key, value in data.items():
        print(f"{key:15s} : {value}")

# -----------------------------------------------------------------------------
if __name__  == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Generate parameters file using JSON format')
    parser.add_argument('-f', '--filename',
                        default='default',
                        type=str, nargs='?',
                        help='Parameters file (json format)')

    parser.add_argument('-po', '--parameters_old',
                        action='store_true',
                        default=False,
                        help='Set this to generate default parameters file instead.')

    parser.add_argument('-o', '--output-directory',
                        nargs='?',
                        type=str,
                        default='configs',
                        help='Output directory')

    args = parser.parse_args()

    # Delete custom connectivity
    if args.parameters_old:
        del _data["connectivity"]["inter_custom"]

    # Make directories
    basedir = args.output_directory
    if not os.path.isdir(basedir):
        print('[+] Creating directory', basedir)
        os.makedirs(basedir)

    # dir_EC = os.path.join(basedir, 'opt_noise_EC')
    # dir_DG = os.path.join(basedir, 'opt_noise_DG')
    # dir_CA3 = os.path.join(basedir, 'opt_noise_CA3')
    # dir_CA1 = os.path.join(basedir, 'opt_noise_CA1')

    vmin,vmax = 30,45
    noise_vals = np.arange(vmin, vmax)

    areas = ['EC', 'DG', 'CA3', 'CA1']
    # for area in areas:
    #     currdir = os.path.join(basedir, 'opt_noise_'+area)
    #     if not os.path.isdir(currdir):
    #         print('[+] Creating directory', currdir)
    #         os.makedirs(currdir)
    #
    #     # Fix the parameters
    #     cnt=0
    #     for val in noise_vals:
    #         _data["areas"][area]["E"]["noise"] = np.around(val*100.e-06, 6)
    #         _data["areas"][area]["I"]["noise"] = np.around(val*10.e-06, 6)
    #
    #         # Define the filename
    #         # filename = "./{0}/{1}_{2:02d}.json".format(args.output_directory, args.filename, cnt)
    #         filename = os.path.join(currdir, (args.filename+'{0:02d}.json').format(cnt))
    #         print('Saving file "{0}"'.format(filename))
    #         save(filename, _data)
    #         cnt += 1
    #
    #         _data["areas"][area]["E"]["noise"] = 0.
    #         _data["areas"][area]["I"]["noise"] = 0.

    vmin,vmax,vstep = 0., 55., 5.
    vals = np.arange(vmin, vmax, vstep)
    # vals = np.linspace(vmin, vmax, 17)
    # stim_times = [5.6050, 5.6153,5.6257, 5.6360, 5.6464, 5.6567, 5.6671, 5.6775, 5.6879, 5.7557, 5.7734, 5.7856, 5.7966, 5.8072, 5.8176, 5.8280, 5.8383]
    # stim_times = np.loadtxt('stim_times.txt')
    # stim_amplitude = 20.
    cnt = 0
    stim_t_off = 1.5
    for val in vals:
        _data["Kuramoto"]["kN"] = np.around(val, 2)
        # _data["fixed_input"]["high"] = np.around(val,2)
        # _data["connectivity"]["inter_custom"]["EC"]["E"][2] = [np.around(val, 2)]*2
        # _data["connectivity"]["inter_custom"]["DG"]["E"][2] = [np.around(val, 2)]*2
        # _data["connectivity"]["inter_custom"]["EC"]["E"][3] = [np.around(val, 2)]*2
        # _data["connectivity"]["inter_custom"]["CA3"]["E"][3] = [np.around(val, 2)]*2
        # _data["connectivity"]["inter_custom"]["CA1"]["E"][0] = [np.around(val, 2)]*2

        # _data["areas"]["EC"]["E"]["noise"] = np.around(val, 6)
        # _data["areas"]["EC"]["I"]["noise"] = np.around(val, 6)
        # _data["areas"]["DG"]["E"]["noise"] = np.around(val, 6)
        # _data["areas"]["DG"]["I"]["noise"] = np.around(val, 6)
        # _data["areas"]["CA3"]["E"]["noise"] = np.around(val, 6)
        # _data["areas"]["CA3"]["I"]["noise"] = np.around(val, 6)
        # _data["areas"]["CA1"]["E"]["noise"] = np.around(val, 6)
        # _data["areas"]["CA1"]["I"]["noise"] = np.around(val, 6)

        # _data["stimulation"]["onset"] = np.round(stim_t_off+val,3)
        # _data["stimulation"]["I"] = [np.round(val,1)]
        # _data["stimulation"]["stim_freq"] = np.around(val, 1).astype(int).item()
        # _data["stimulation"]["nr_of_trains"] = np.rint(5*val).astype(int).item()

        # Used for PRC calc
        # _data["stimulation"]["I"] = [np.round(stim_amplitude,1)]
        # _data["stimulation"]["onset"] = np.round(val,4)

        # Define the filename
        filename = os.path.join(basedir, (args.filename+'{0:02d}.json').format(cnt))
        print('Saving file "{0}"'.format(filename))
        save(filename, _data)
        cnt += 1

    # _data["connectivity"]["inter_custom"]["EC"]["E"][1] = np.around([a, a], decimals=1).tolist()
    # _data["connectivity"]["inter_custom"]["EC"]["E"][2] = np.around([b, b], decimals=1).tolist()
    # _data["connectivity"]["inter_custom"]["EC"]["E"][3] = np.around([c, c], decimals=1).tolist()
    # _data["connectivity"]["inter_custom"]["DG"]["E"][2] = np.around([b, b], decimals=1).tolist()
    # _data["connectivity"]["inter_custom"]["CA3"]["E"][3] = np.around([c, c], decimals=1).tolist()
    # _data["connectivity"]["inter_custom"]["CA1"]["E"][0] = np.around([d, d], decimals=1).tolist()
