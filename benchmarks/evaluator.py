def get_component_value(vec, comp_idx):
    """
    vec: data vector
    comp_idx: int or iterable of ints
    """
    if isinstance(comp_idx, (list, tuple)):
        return sum(vec[i] for i in comp_idx if i < len(vec))
    else:
        return vec[comp_idx] if comp_idx < len(vec) else 0.0


class Evaluator:
    def __init__(self, experiments, base_setting, compare_setting, components):
        self.experiments = experiments
        self.base_setting = base_setting     # "s0"
        self.compare_setting = compare_setting  # "s1"
        self.components = components

    def admissible(self, l, p):
        return p <= 4 and (3*l - p) >= 17
#        return p <= 5 and (3*l - p) >= 16

    def evaluate(self, summary_data, output_file="improvement_table.tex"):
        """
        Produces a LaTeX table:
        - rows: geometries / experiments
        - columns: components
        Uses exactly three (l,p) points:
            (8,5), (9,8), (10,11)
        Displays averaged improvement in percent, rounded to integer.
        """

        geometries = sorted({k[0] for k in summary_data.keys()})

        # fixed admissible points
        points = [(8, 5), (8, 6), (8, 7), (9, 8), (9, 9), (9, 10), (10, 11)]

        with open(output_file, "w") as f:
            # ---------- LaTeX header ----------
            cols = "l" + "c" * len(self.components)
            f.write("\\begin{tabular}{" + cols + "}\n")
            f.write("\\hline\n")

            header = ["Geometry"] + [comp[1].capitalize() for comp in self.components]
            f.write(" & ".join(header) + " \\\\\n")
            f.write("\\hline\n")

            # ---------- Table body ----------
            for e_folder in geometries:
                row = [self.experiments[e_folder]]
                row_adj = [self.experiments[e_folder] + ",adj"]

                for comp in self.components:
                    comp_idx = comp[0]

                    improvements = []
                    improvements_adjusted = []

                    for l, p in points:
                        try:
                            base = summary_data[(e_folder, self.base_setting, l, p)]
                            comp = summary_data[(e_folder, self.compare_setting, l, p)]
                        except KeyError:
                            continue

                        val1 = get_component_value(base,comp_idx)
                        val0 = get_component_value(comp,comp_idx)

                        if val1 == 0:
                            continue

                        imp = (val1 - val0) / val1
                        improvements.append(imp)

                        # adjusted improvement
                        try:
                            num_idx = 7
                            num1 = base[num_idx]
                            num0 = comp[num_idx]
                        except KeyError:
                            continue

                        if num1 != num0 and num1 != 0:
                            adj_val1 = val1 * num0 / num1
                            adj_imp = (adj_val1 - val0) / adj_val1
                            improvements_adjusted.append(adj_imp)

                    # ---------- averaging ----------
                    if improvements:
                        avg = round(sum(improvements) / len(improvements), 2)
                        row.append(f"{avg:.2f}")
                    else:
                        row.append("--")

                    if improvements_adjusted:
                        avg_adj = round(sum(improvements_adjusted) / len(improvements_adjusted), 2)
                        row_adj.append(f"{avg_adj:.2f}")

                f.write(" & ".join(row) + " \\\\\n")
                if len(row_adj) > 1:
                    f.write(" & ".join(row_adj) + " \\\\\n")

            # ---------- Footer ----------
            f.write("\\hline\n")
            f.write("\\end{tabular}\n")
