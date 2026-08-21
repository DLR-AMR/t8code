# ------------------------------------------------------------
# plotter.py — version with grouped-per-component figures
# ------------------------------------------------------------
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

def get_component_value(vec, comp_idx):
    """
    vec: data vector
    comp_idx: int or iterable of ints
    """
    if isinstance(comp_idx, (list, tuple)):
        return sum(vec[i] for i in comp_idx if i < len(vec))
    else:
        return vec[comp_idx] if comp_idx < len(vec) else 0.0


class Plotter:
    def __init__(self, experiments, settings, components, colorscheme, root="."):
        self.experiments = experiments
        self.settings = settings
        self.components = components
        self.colorscheme = colorscheme
        self.root = root

    # --------------------------------------------------------
    # MASTER PLOT FUNCTION
    # --------------------------------------------------------
    def plot(self, summary, components, formats, out_root):
#         geometries = sorted(self.experiments.keys())
        geometries = self.experiments.keys()
        n_geoms = len(geometries)

        # 2-column layout
        ncols = 2
        nrows = math.ceil(n_geoms / ncols)

        for component in components:
            comp_idx, comp_name, factor, weak_scaling, strong_scaling, scaling_lower_bound = component
            fig, axes = plt.subplots(
                nrows=nrows, ncols=ncols,
                figsize=(10, 4*nrows),
                squeeze=True,
                sharex=True,
                sharey=True
            )

            # prepare color map once per component
            levels = sorted({l for (e,s,l,p) in summary.keys()})
            color_map = {l: self.colorscheme(i) for i,l in enumerate(levels)}

            # fill subplots
            for k, e_folder in enumerate(geometries):
                ax = axes[k // ncols][k % ncols]

                legend_items = self._plot_single_geometry_into_axis(
                    ax, summary, e_folder, component, color_map)

                ax.set_title(self.experiments.get(e_folder, e_folder))

            # hide empty axes if odd number of geometries
            for k in range(n_geoms, nrows*ncols):
                axes[k//ncols][k%ncols].axis("off")

            # ----------------------------------------------------
            # SINGLE SHARED LEGEND BELOW ALL SUBPLOTS
            # ----------------------------------------------------
            handles, labels = legend_items
            fig.legend(
                handles, labels,
                loc="lower center",
                ncol=(len(handles)+1)//2,
                fontsize=9,
                markerscale=0.9,
                bbox_to_anchor=(0.5, -0.03)
            )

            fig.suptitle(comp_name.capitalize(), fontsize=14)
            fig.tight_layout(rect=[0, 0.05, 1, 0.97])

            # save only ONE figure per component
            self._save_grouped(fig, comp_name, formats, out_root)
            plt.close(fig)

    # --------------------------------------------------------
    # Helper: save grouped figure
    # --------------------------------------------------------
    def _save_grouped(self, fig, comp_name, formats, out_root):
        base = comp_name
        for ext in formats:
            out_dir = os.path.join(out_root, ext)
            os.makedirs(out_dir, exist_ok=True)
            fig.savefig(os.path.join(out_dir, base + "." + ext), bbox_inches="tight")

    # --------------------------------------------------------
    # Plot one geometry into an axis
    # --------------------------------------------------------
    def _plot_single_geometry_into_axis(
        self, ax, summary, e_folder, component, color_map
    ):
        comp_idx, comp_name, factor, weak_scaling, strong_scaling, scaling_lower_bound = component

        # axes, ticks, limits
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")

        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        formatter.set_useOffset(False)
        ax.xaxis.set_major_formatter(formatter)

        ax.set_xlabel("Number of processes")
        ax.set_ylabel("Time")

        ymin, ymax = 1e-4 / factor, 30/factor
        ax.set_ylim(ymin, ymax)

        ax.grid(True, which="major", linestyle="--", alpha=0.3)
        ax.grid(True, which="minor", linestyle=":", alpha=0.3)
        ax.minorticks_on()

        # level data
        keys_sorted = sorted(summary.keys())
        for (e, s, l, p) in keys_sorted:
            if e != e_folder:
                continue

            vec = summary[(e,s,l,p)]

            y = get_component_value(vec, comp_idx)

            x = 2**p
            c = color_map.get(l)

            if s == "s0":
                ax.plot(x, y, linestyle="none", marker="x", markersize=8, color=c)
            else:
                ax.scatter(x, y, marker="s", s=64, facecolors="none", edgecolors=c)

        legend_levels = [
            (Line2D([0],[0], color=color_map[l], marker="o",
                    linestyle="None", markersize=6), f"Level {l}")
            for l in sorted(color_map.keys())
        ]
        legend_shapes = [
            (Line2D([0],[0], color="black", marker="x", linestyle="None", markersize=8), "New"),
            (Line2D([0],[0],
                    color="black",
                    marker="s",
                    linestyle="None",
                    markersize=8,
                    markerfacecolor="none")
                , "Old")
            ]

        items =legend_levels + legend_shapes


        # strong scaling
        any_strong = False
        any_superstrong = False
        for s,l,p, superstrong in strong_scaling:       
            any_strong = True
            y0 = get_component_value(summary[(e_folder,s,l,p)], comp_idx)

            if superstrong:
                C = y0 * (2**p) / (3*l-p)
                xs = np.linspace(2**p, 2**15, 300)
                ys = np.array([C/x*(3*l-math.log2(x)) for x in xs])
                ax.plot(xs, ys, color=color_map[l], linestyle=(0, (1, 8)), linewidth=1.2)
                any_superstrong = True

            C = y0 * (2**p)
            xs = np.linspace(2**p, 2**15, 300)
            ax.plot(xs, C/xs, color=color_map[l], linestyle=":", linewidth=1.2)
        
        if any_strong:
            items.append((Line2D([0],[0], color=color_map[l], linestyle=":", linewidth=1.2),
                                    f"Strong Scaling"))
        if any_superstrong:
            items.append((Line2D([0],[0], color="black", linestyle=(0, (1, 8)), linewidth=1.2),
                                              f"Super Scaling"))



        # weak scaling
        any_weak = False
        for s,l,p in weak_scaling:       
            any_weak = True
            y0 = get_component_value(summary[(e_folder,s,l,p)], comp_idx)
            xs = np.linspace(2**p, 2**15, 300)
            ax.plot(xs, [y0]*len(xs), color="gray", linestyle="--", linewidth=1.2)
        if any_weak:
            items.append((Line2D([0],[0], color="gray", linestyle="--", linewidth=1.2), "Weak Scaling"))


        # scaling bound
        if scaling_lower_bound:
            pvals = sorted(scaling_lower_bound[e].keys())
            xs = np.array([2**p for p in pvals])
            ys_min = np.array([scaling_lower_bound[e][p] for p in pvals])
            ax.plot(xs, ys_min, color="black", linewidth=2.0)
            items.append                ((Line2D([0],[0], color="black", linewidth=2), "Scaling Bound"))

        handles = [h for h,_ in items]
        labels  = [t for _,t in items]
        return handles, labels
