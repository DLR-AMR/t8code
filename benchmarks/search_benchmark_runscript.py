from pathlib import Path
import subprocess

eclasses = [4, 5, 6, 7]
schemes = [0, 1]
max_proc_log = 14

max_elem_log=20
min_elem_log=5
repititions = 20

maxlevel=18
dim=3

#CARO has 128 =2^^7 procs per node
log_procs_per_node=7

out_directory = "out"

for eclass in eclasses:
    for scheme in schemes:
        for proc_log in range(max_proc_log+1):
                n_nodes = max(1, int(2**(proc_log- log_procs_per_node)))
                string="#!/bin/bash\n"
                string += "#SBATCH --time=00:15:00\n"
                string += f"#SBATCH -N {n_nodes}\n"
                string += f"#SBATCH --ntasks={2**proc_log}\n"
                for level in range(maxlevel+1,-1,-1):
                    if (max_elem_log >= dim*level - proc_log >= min_elem_log):
                        string += f"echo srun ./t8_time_search -s{scheme} -e{eclass} -l{level} -r{repititions}\n"
                        string += f"srun -o {out_directory}/e{eclass}/s{scheme}/p{proc_log}_l{level}_%j ./t8_time_search -s{scheme} -e{eclass} -l{level} -r{repititions}\n"

                file = Path(f"{out_directory}/e{eclass}/s{scheme}/scripts/p{proc_log}.sh")
                file.parent.mkdir(exist_ok=True, parents=True)
                file.write_text(string)

for eclass in eclasses:
    for scheme in schemes:
        for proc_log in range(max_proc_log+1):
            subprocess.call(f"sbatch {out_directory}/e{eclass}/s{scheme}/scripts/p{proc_log}.sh", shell=True)
