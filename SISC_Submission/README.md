# Instructions for benchmark execution of SISC Paper 
## "A MORTON-FREUDENTHAL-KUHN SPACE-FILLING CURVE FOR HYBRID ADAPTIVE MESH REFINEMENT"


First checkout and build this t8code branch:
```
git clone --branch SISC_benchmark_MFK --depth 1 git@github.com:DLR-AMR/t8code.git
cd t8code
git submodule update --init
mkdir build && cd build
cmake ..
make -j
cd benchmarks
```
Now you can execute single instances of the t8_time_search benchmark with
```
./t8_time_search -s<scheme> -l<level> -e<element> -r<repitition>
```
where the options have the following meaning:

```
[t8]    -s | --scheme    <INT>       Option to choose the scheme, 0: standalone scheme, 1: default scheme
[t8]    -l | --level     <INT>       The level of the forest.
[t8]    -e | --elements  <INT>       Specify the type of elements to use.
					2 - quadrilateral
					3 - triangle
					4 - hexahedron (default)
					5 - tetrahedron
					6 - prism
					7 - pyramid
[t8]    -r | --repetitions <INT>     The number of repeats of the new-search-adapt-partition cycle.
```

You can also submit batch jobs on slurm with:

```
python3 search_benchmark_runscript.py
```
This creates an appropriate `.batch` file for each configuration and submits it. You might want to start this script on a smaller number of maximal ranks by reducing `log_procs_per_node` inside the script, or reduce the maximal effort per rank by reducing `max_elem_log`.

If you submitted the batch files via the runscript, the result files should be in correct directory structure so that you can use 

```
python3 search_benchmark_visscript.py
```
to create the visualizations and runtime improvement table from the paper.
