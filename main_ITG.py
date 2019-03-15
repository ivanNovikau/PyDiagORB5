import Mix as mix
import ITG_gamma as itg
import numpy as np

# Main Script: to plot non-zonal components,
# calculate growth rates


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(itg)


def prepare_data(dd, path_to_folder):
    dd = itg.read_data(dd, path_to_folder)
    return dd


path_sim = 'd:/Work-Projects/MyProgs/ORB_data/' \
                     'turbulence/AUG20787/adiab/linear/n60/'
                     
                     
path_sim_draco = '/u/ivannovi/turbulence/AUG20787/adiab/linear/n60'

# Main principles:
# -> main script which import necessary modules
# -> some modules can have functions with very close functionality, but this fact will
# allow to keep the diagnostic more flexible
# -> i can keep a local folder (on a given machine) with some global variables:
# paths to created projects, their names, list of variables, description of a project and
# its variables, project configuration
# -> it should be possible to use main scripts without creating any project, and within
# a particular project
# -> Using main script without a project, one can save data to temporal files
# (to avoid overloading of the operational memory)
# -> With a project one can define its configuration, save data to permanent files to load
# these data afterwards

# Structure
# Project -> Main Script -> modules
# or
# Main Script -> modules

# Questions:
# How to set project configuration, that would be specific for a particular project without
# making complex the whole project?


