import Mix as mix
import ITG_gamma as itg
import zf_gam as zf
import numpy as np

# Main Script: to plot non-zonal components,
# calculate growth rates


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(itg)
    mix.reload_module(zf)


path_n60 = 'd:/Work-Projects/MyProgs/ORB_data/' \
                     'turbulence/AUG20787/adiab/linear/n60/'
path_sim_draco = '/u/ivannovi/turbulence/AUG20787/adiab/linear/n60'

dd = {}
dd['path'] = path_n60

                     




