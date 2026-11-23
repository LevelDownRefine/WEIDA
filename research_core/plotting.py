"""Common descriptive distribution and paired-change figures for every backend."""
from pathlib import Path
import numpy as np


def plot_study(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .study import cohort_ids

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    ids = cohort_ids(result.scans)
    indices = {g: sorted([i for i, m in enumerate(result.scans) if m['condition'] == g],
                         key=lambda i: ids[i]) for g in ('HC', 'OFF', 'ON')}
    colors = ['#7c899b', '#d58b40', '#5b9d9a']
    algorithm = result.algorithm
    for family_index, (family, values) in enumerate(algorithm.features.items(), 1):
        for start in range(0, values.shape[1], 16):
            columns = list(range(start, min(start+16, values.shape[1])))
            ncols = min(4, len(columns))
            nrows = (len(columns) + ncols-1) // ncols
            for paired in (False, True):
                fig, axes = plt.subplots(nrows, ncols, squeeze=False, figsize=(3.3*ncols, 2.8*nrows))
                for column, ax in zip(columns, axes.flat):
                    if paired:
                        for a, b in zip(values[indices['OFF'], column], values[indices['ON'], column]):
                            if np.isfinite([a, b]).all():
                                ax.plot([0, 1], [a, b], color=colors[2], alpha=.5, lw=.8)
                        ax.set_xticks([0, 1], ['OFF', 'ON'])
                    else:
                        for position, group in enumerate(('HC', 'OFF', 'ON')):
                            y = values[indices[group], column]
                            y = y[np.isfinite(y)]
                            if len(y) > 1 and np.ptp(y) > 0:
                                ax.violinplot(y, positions=[position], showmeans=True, showextrema=False)
                            ax.scatter(np.full(len(y), position), y, color=colors[position], s=9)
                        ax.set_xticks([0, 1, 2], ['HC', 'OFF', 'ON'])
                    ax.set_title(algorithm.names[family][column])
                    ax.set_ylabel(algorithm.units[family])
                for ax in list(axes.flat)[len(columns):]:
                    ax.set_visible(False)
                fig.suptitle(f'{algorithm.method} | {family} | descriptive')
                fig.tight_layout()
                kind = 'paired_changes' if paired else 'group_distributions'
                for extension in ('png', 'pdf'):
                    fig.savefig(output / f'{kind}_family{family_index}_{start//16+1}.{extension}', dpi=160)
                plt.close(fig)
