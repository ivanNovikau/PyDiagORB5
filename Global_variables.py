from IPython.core.getipython import get_ipython


def reload():
    return


# ---------------------------------------------------------------------------
MIN_N_PEAKS = 3
COEF_ERR = 1.96          # alpha-quantile of the standard normal distribution
CONFIDENCE_PERC = 0.95   # (2*alpha - 1)-confidence interval, that corresponds to the alpha quantile

# ---------------------------------------------------------------------------
NONE_FILTER = {'sel_filt': None}

# ---------------------------------------------------------------------------
if 'Terminal' in get_ipython().__class__.__name__:
    FLAG_LATEX = True
else:
    FLAG_LATEX = False

if FLAG_LATEX:
    FIG_SIZE_W = 15
    FIG_SIZE_H = 9.5
    LEG_SCALE = 1.0
    FONT_SIZE = 28
    FLAG_LATEX = True
    FONT_SIZE_LABELS = FONT_SIZE * 1.8
    FONT_SIZE_TICKS = FONT_SIZE * 1.5
    FONT_SIZE_ORDER = FONT_SIZE * 1.2
    FONT_SIZE_TITLE = FONT_SIZE * 1.3

    # # configuration #2.1 (used for q-plot in AAPPS paper)
    # FIG_SIZE_W = 9.5
    # FIG_SIZE_H = 9.5
    # LEG_SCALE = 1.0
    # FONT_SIZE = 28
    # FLAG_LATEX = True
    # FONT_SIZE_LABELS = FONT_SIZE * 1.8
    # FONT_SIZE_TICKS = FONT_SIZE * 1.2
    # FONT_SIZE_ORDER = FONT_SIZE * 1.2
    # FONT_SIZE_TITLE = FONT_SIZE * 1.3

    # # configuration #3 (used for B plot in AAPPS paper)
    # FIG_SIZE_W = 9.5
    # FIG_SIZE_H = 9.5
    # LEG_SCALE = 0.82
    # FONT_SIZE = 28
    # FLAG_LATEX = True
    # FONT_SIZE_LABELS = FONT_SIZE * 2
    # FONT_SIZE_TICKS = FONT_SIZE * 1.2
    # FONT_SIZE_ORDER = FONT_SIZE * 2
    # FONT_SIZE_TITLE = FONT_SIZE * 1.5
else:
    FIG_SIZE_W = 10
    FIG_SIZE_H = 6
    LEG_SCALE = 0.5
    FONT_SIZE = 22
    FLAG_LATEX = False
    FONT_SIZE_LABELS = FONT_SIZE * 0.6
    FONT_SIZE_TICKS  = FONT_SIZE * 0.6
    FONT_SIZE_ORDER  = FONT_SIZE * 0.6
    FONT_SIZE_TITLE  = FONT_SIZE * 0.6

# ---------------------------------------------------------------------------
# --- DEFAULT VARIABLE DEFINITIONS ---
# ---------------------------------------------------------------------------
def_erbar_ts = {
    'type':             'zonal',
    'variable':         'er',
    'plane':            'ts',
    'avr_operation':    'point-s',
    'avr_domain':       0.5,
}
def_safety_factor = {
    'type': 'equ-profile',
    'variable': 'q',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}
def_Teq_deuterium = {
    'type': 'equ-profile',
    'variable': 'T-equ',
    'species_name': 'deuterium',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}


def create_signal(default_signal, dd):
    res_signal = dict(default_signal)
    res_signal.update({
        'dd': dd
    })
    return [res_signal]