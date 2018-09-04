import matplotlib as mpl
import matplotlib.pyplot as plt

output = 'pdf'


def figsize(scale):
    fig_width_pt = 469.75499
    inches_per_pt = 1.0 / 72.27
    # golden_mean = (np.sqrt(5.0) - 1.0) / 2.0
    golden_mean = 0.7
    fig_width = fig_width_pt * inches_per_pt * scale
    fig_height = fig_width * golden_mean
    fig_size = [fig_width, fig_height]
    return fig_size


def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax


def savepgf(filename):
    plt.tight_layout()
    plt.savefig('output/{}.{}'.format(filename, output))


pgf_with_latex = {
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,
    "font.family": "Latin Modern Roman",
    "font.serif": [],
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(1),
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}", r"\usepackage[T1]{fontenc}",
        r"\usepackage{lmodern}", r"\usepackage{textcomp, gensymb}"
    ]
}

mpl.rcParams.update(pgf_with_latex)
