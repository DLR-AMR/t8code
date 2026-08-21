if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from data_manager import DataManager
    from grouped_plot import Plotter
    from evaluator import Evaluator

    experiments = {"e4": "hexahedron", "e5": "tetrahedron", "e6": "prism", "e7": "pyramid"}
    settings = {"s0": "standalone", "s1": "default"}

    dm = DataManager(experiments, settings)

    summary = dm.load("Summary")

    scaling_limit_idx = 9
    scaling_lower = dm.extract_scaling_boundary(summary, scaling_limit_idx, min)

    scaling_lower_total = {k: {p: v*3 for p,v in d.items()} for k,d in scaling_lower.items()}

    components = [
        [0, "total", 1, [("s0",5,0),("s0",7,1)], [("s0",5,0, False),("s0",7,1,False)], scaling_lower_total],
        [[1,3,4], "AMR", 1, [("s0",8,4)], [("s0",8,4,False)], scaling_lower_total],
        [1, "new",10, [("s0",7,1)],[("s0",7,1, False)],scaling_lower],
        [3, "refine", 10, [("s0",7,1)],[("s0",7,1, False)],scaling_lower],
        [4, "partition",10, [("s0",8,4)], [("s0",8,4, False)], scaling_lower],
        [2, "search",1,[("s0",5,0),("s0",7,1)], [("s0",5,0, True),("s0",7,1,True)],None]
    ]

    colorscheme = lambda i: plt.cm.tab20(i % 20)

    plotter = Plotter(experiments, settings, components, colorscheme)
    plotter.plot(summary, components, ["pdf"], "out")

    evaluator = Evaluator(experiments, "s1", "s0", components)
    evaluator.evaluate(summary)
