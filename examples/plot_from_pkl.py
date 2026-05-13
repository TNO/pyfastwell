"""
Standalone plotting script that reads model_io.pkl and regenerates:
  - stoch_expt  (expectation / exceedance) plots
  - stoch_cross (cross-plot) plots, with dots coloured by a chosen variable

Colouring options
-----------------
  --colorby k        : 5 equally-spaced k-value classes   (original behaviour)
  --colorby power    : fixed target-power classes at 4, 6, 8, 10, 12, 14 MW
                       (values below 4 -> class 0, above 14 -> last class)

No model re-run required.
"""

import os
import pickle
import argparse

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ── output-array index constants (mirror fastmodelBLT01_core.py) ────────────
ITIME      = 0
IPOWER     = 1
ICOP       = 2
INPV       = 3
ILCOH      = 4
ITEMPPRD   = 5
IFLOWRATE  = 6

RESNAMES  = ['POWER [MW]', 'COP [-]', 'NPV [mln EUR]',
             'LCOE [EUR MWh-1]', 'TEMPPRD [C]', 'FLOWRATE [m3 h-1]']
INDICES   = [IPOWER, ICOP, INPV, ILCOH, ITEMPPRD, IFLOWRATE]
PARAM_UNITS = ['[mDarcy]', '[m]']

# Fixed power class edges requested: 4, 6, 8, 10, 12, 14  (step=2)
POWER_EDGES = list(range(4, 15, 2))   # [4, 6, 8, 10, 12, 14]


# ── helpers ─────────────────────────────────────────────────────────────────

def ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path)


def _make_colours(n_classes):
    """Return a list of RGBA colours for n_classes using the viridis colormap."""
    # Use pyplot.get_cmap(name, lut) for broad Matplotlib version compatibility.
    cmap = plt.get_cmap('viridis', n_classes)
    return [cmap(i / max(n_classes - 1, 1)) for i in range(n_classes)]


def k_classes(k_values, n_classes=5):
    """
    Assign each realisation to one of *n_classes* equally-spaced k-value bins.
    Returns (class_indices 1-based, bin_edges, colours, labels).
    """
    edges   = np.linspace(k_values.min(), k_values.max(), n_classes + 1)
    indices = np.digitize(k_values, edges, right=True)
    indices = np.clip(indices, 1, n_classes)
    colours = _make_colours(n_classes)
    labels  = [f'k {edges[i]:.0f}-{edges[i+1]:.0f} mDarcy' for i in range(n_classes)]
    return indices, n_classes, colours, labels


def power_classes(power_values, edges=None):
    """
    Assign each realisation to a power class defined by *edges*.

    The edges list [e0, e1, ..., eN] defines N+1 bins:
        class 1  : power < e0
        class 2  : e0 <= power < e1
        ...
        class N+1: power >= eN

    Returns (class_indices 1-based, n_classes, colours, labels).
    """
    if edges is None:
        edges = POWER_EDGES

    n_classes = len(edges) + 1          # one extra bin on each outer side
    colours   = _make_colours(n_classes)

    # Build labels
    labels = [f'< {edges[0]} MW']
    for i in range(len(edges) - 1):
        labels.append(f'{edges[i]}-{edges[i+1]} MW')
    labels.append(f'>= {edges[-1]} MW')

    # Assign class (1-based)
    indices = np.ones(len(power_values), dtype=int)
    for i, edge in enumerate(edges):
        indices[power_values >= edge] = i + 2   # classes 2..n_classes
    indices = np.clip(indices, 1, n_classes)

    return indices, n_classes, colours, labels



# ── expectation plot ─────────────────────────────────────────────────────────

def expectation_plot(results, name, pvals=None, filename_noext=None):
    """Reproduce Stochasticmodel.expectation_plot()."""
    tofile     = f'{filename_noext}_expt_{name}.png' if filename_noext else None
    tofile_csv = f'{filename_noext}_expt_{name}.csv' if filename_noext else None

    if tofile:
        ensure_dir(os.path.dirname(tofile))

    sorted_samples = np.sort(results)[::-1]          # descending → exceedance
    cdf = np.arange(1, len(sorted_samples) + 1) / len(sorted_samples) * 100

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(sorted_samples, cdf, linestyle='-')
    ax.set_xlabel(name)
    ax.set_ylabel('Cumulative Density [%]')
    ax.set_ylim(0, 100)
    ax.set_title(f'Expectation plot of {name}')
    ax.grid(True)

    if pvals is not None:
        for pval in pvals:
            v = np.interp(pval, cdf, sorted_samples)
            ax.plot(v, pval, 'ro')
            label = f'p{pval} = {v:.2f}'
            ax.annotate(label, (v, pval), textcoords='offset points',
                        xytext=(5, 5), ha='center')

    if tofile:
        fig.savefig(tofile, dpi=150, bbox_inches='tight')
        df = pd.DataFrame({name: sorted_samples, 'cdf': cdf})
        df.to_csv(tofile_csv, index=False)
        print(f'  saved: {tofile}')
    else:
        plt.show()
    plt.close(fig)


# ── cross-plot with class colouring ─────────────────────────────────────────

def cross_plot_colored(inputval, inputname, results, name,
                       class_idx, n_classes, colours, labels,
                       color_varname, color_varvals,
                       filename_noext=None):
    """
    Cross-plot of inputval vs results.  Dots are coloured by the supplied
    class assignments (either k-based or power-based).

    Parameters
    ----------
    inputval      : x-axis values (sampled parameter)
    inputname     : x-axis label
    results       : y-axis values (model output)
    name          : y-axis label
    class_idx     : 1-based integer class per realisation
    n_classes     : total number of classes
    colours       : list of RGBA colours, length n_classes
    labels        : list of class legend labels, length n_classes
    color_varname : name of the colouring variable (for title / CSV)
    color_varvals : raw values of the colouring variable (for CSV)
    filename_noext: path without extension; None -> show interactively
    """
    tofile     = f'{filename_noext}_cross_{name}_{inputname}.png' if filename_noext else None
    tofile_csv = f'{filename_noext}_cross_{name}_{inputname}.csv' if filename_noext else None

    if tofile:
        ensure_dir(os.path.dirname(tofile))

    fig, ax = plt.subplots(figsize=(10, 6))

    for cls in range(1, n_classes + 1):
        mask = class_idx == cls
        if mask.any():
            ax.scatter(inputval[mask], results[mask],
                       color=colours[cls - 1], label=labels[cls - 1],
                       alpha=0.75, edgecolors='none', s=25)

    ax.set_xlabel(inputname)
    ax.set_ylabel(name)
    ax.set_title(f'{name} vs {inputname}  (coloured by {color_varname})')
    ax.legend(title=color_varname, fontsize=8, title_fontsize=8,
               loc='best', framealpha=0.7)
    ax.grid(True)

    if tofile:
        fig.savefig(tofile, dpi=150, bbox_inches='tight')
        df = pd.DataFrame({inputname: inputval, name: results,
                           color_varname: color_varvals,
                           'color_class': class_idx})
        df.to_csv(tofile_csv, index=False)
        print(f'  saved: {tofile}')
    else:
        plt.show()
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Re-generate stoch_expt and stoch_cross plots from model_io.pkl')
    parser.add_argument('pkl', nargs='?',
                        default='output/BLT01_100mDarcy/model_io.pkl',
                        help='Path to model_io.pkl  (default: output/BLT01_100mDarcy/model_io.pkl)')
    parser.add_argument('--outdir', '-o', default=None,
                        help='Output directory for plots (default: same folder as pkl)')
    parser.add_argument('--colorby', '-c', default='k',
                        choices=['k', 'power'],
                        help='Variable used to colour cross-plot dots: '
                             '"k" = 5 equal k classes (default); '
                             '"power" = fixed power classes 4,6,8,10,12,14 MW')
    parser.add_argument('--nclasses', '-n', type=int, default=5,
                        help='Number of k-value colour classes when --colorby k (default: 5)')
    args = parser.parse_args()

    pkl_path       = args.pkl
    outdir         = args.outdir or os.path.dirname(os.path.abspath(pkl_path))
    filename_noext = os.path.join(outdir, 'stoch')

    print(f'Loading {pkl_path} ...')
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)

    model = data['model']
    M     = model.stoch.M       # list of param arrays, len = nsamples
    fwi   = model.stoch.fwi     # shape (nsamples, nresults, ntimes)

    M_arr       = np.array(M)                            # (nsamples, nfree_params)
    free_params = model.stoch.argdict['free_param']      # e.g. ['k', 'Lh']

    # ── build colour classes ────────────────────────────────────────────────
    if args.colorby == 'power':
        power_vals = fwi[:, IPOWER, ITIME]
        class_idx, n_classes, colours, labels = power_classes(power_vals, POWER_EDGES)
        color_varname = 'POWER [MW]'
        color_varvals = power_vals
        print(f'Colouring by POWER with edges {POWER_EDGES} MW  '
              f'({n_classes} classes)')
    else:  # 'k'
        try:
            k_idx = free_params.index('k')
        except ValueError:
            print("WARNING: 'k' not found in free parameters; using first parameter.")
            k_idx = 0
        k_vals = M_arr[:, k_idx]
        class_idx, n_classes, colours, labels = k_classes(k_vals, args.nclasses)
        color_varname = 'k [mDarcy]'
        color_varvals = k_vals
        print(f'Colouring by k  ({n_classes} equal classes)')

    pvals = [10, 30, 50, 70, 90]

    print('\nGenerating expectation plots ...')
    for ires, resname in enumerate(RESNAMES):
        index = INDICES[ires]
        expectation_plot(fwi[:, index, ITIME], resname,
                         pvals=pvals, filename_noext=filename_noext)

    print('\nGenerating cross-plots ...')
    # Build axis labels for free parameters
    unit_map = {'k': '[mDarcy]', 'Lh': '[m]'}
    for ires, resname in enumerate(RESNAMES):
        index = INDICES[ires]
        for i, param in enumerate(free_params):
            unit = unit_map.get(param, '')
            xlabel = f'{param} {unit}'.strip()
            cross_plot_colored(
                M_arr[:, i], xlabel,
                fwi[:, index, ITIME], resname,
                class_idx, n_classes, colours, labels,
                color_varname, color_varvals,
                filename_noext=filename_noext
            )

    print('\nDone.')


if __name__ == '__main__':
    main()

