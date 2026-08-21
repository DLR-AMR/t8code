import os
import re
import numpy as np
from collections import defaultdict

class DataManager:
    def __init__(self, experiments, settings, root="."):
        self.experiments = experiments      # {"e4": "cube", ...}
        self.settings = settings            # {"s0": "standalone", ...}
        self.root = root

        self.file_regex = re.compile(r"p(?P<p>\d+)_l(?P<l>\d+)_\d+")
        self.line_regex_template = r"{kind}\s*=\s*\[\s*(.*?)\s*\];"

    # ------------------------------------------------------------
    # LOAD
    # ------------------------------------------------------------
    def load(self, kind):
        """
        Returns dict keyed by (e_folder, s_folder, l, p) = avg_vec.
        """
        line_regex = re.compile(self.line_regex_template.format(kind=kind))
        data = {}   # key → avg_vector

        for e_folder in self.experiments.keys():
            e_path = os.path.join(self.root, e_folder)
            if not os.path.isdir(e_path):
                continue

            for s_folder in self.settings.keys():
                s_path = os.path.join(e_path, s_folder)
                if not os.path.isdir(s_path):
                    continue

                for fname in os.listdir(s_path):
                    m = self.file_regex.match(fname)
                    if not m:
                        continue

                    p = int(m.group("p"))
                    l = int(m.group("l"))

                    path = os.path.join(s_path, fname)

                    vectors = []
                    with open(path, "r") as f:
                        for line in f:
                            mm = line_regex.search(line)
                            if mm:
                                vec = np.array(list(map(float, mm.group(1).split())))
                                vectors.append(vec)

                    if not vectors:
                        continue

                    avg = np.mean(vectors, axis=0)
                    key = (e_folder, s_folder, l, p)
                    if not (key in data and avg[0] >= data[key][0]):
                        data[key]=avg

        return data

    # ------------------------------------------------------------
    # Extract min/max per p over all l
    # ------------------------------------------------------------
    def extract_scaling_boundary(self, data, idx, reducer):
        """
        Returns dict: (e_folder → dict p → min/max avg_vec[idx] over all l)
        """
        out = defaultdict(dict)

        for (e, s, l, p), vec in data.items():
            val = vec[idx]
            if p not in out[e]:
                out[e][p] = val
            else:
                out[e][p] = reducer(out[e][p], val)

        return out

    # ------------------------------------------------------------
    # Extract per-level values
    # ------------------------------------------------------------
    def extract_per_level(self, data, e_folder, comp_idx):
        """
        Returns dict: (l → dict p → value)
        """
        per_level = defaultdict(dict)

        for (e, s, l, p), vec in data.items():
            if e == e_folder and comp_idx < len(vec):
                per_level[l][p] = vec[comp_idx]

        return per_level
