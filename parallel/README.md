# motif-based community detection

## Sequential Python implementation

Check the notebook in the sequential_python folder.

## Parallel C++ implementation

The framework is built based on GBBS[https://github.com/ParAlg/gbbs]. The config setup follows GBBS and uses Homegrown scheduler. 
* g++ &gt;= 9.4.0 with support for Cilk Plus
* g++ &gt;= 9.4.0 with pthread support (Homemade Scheduler)

Compile with bazel[https://bazel.build/install]:
```
cd parallel
bazel build //benchmarks/community/wt:Thresholding_main
```
To run
```
bazel run //benchmarks/community/wt:Thresholding_main -- -rounds 1 -community absolute/path/to/groundtruth/community -method METHOD -threshold DELTA -s absolute/path/to/graph/adj
```
Parameters: 
- **rounds** set the number of repetitions.
- **method**, TW, TECTONIC, JACCARD and K3.
- **threshold** the delta value.
- **-s** indicating the input graph is symmetric.

For example, the following command runs TW on DBLP with threshold -15
```
bazel-3.7.0 run //benchmarks/community/wt:Thresholding_main -- -rounds 1 -community path_to_root/large_graphs/dblp.all -method TW -threshold -15.0 -s path_to_root/large_graphs/dblp.adj
```
Available motif-based similarity methods are: TW, TECTONIC, JACCARD and K3.

## Input Formats
-----------
Following GBBS, we use the adjacency graph format that is also supported by the [Problem Based Benchmark
suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html)
and [Ligra](https://github.com/jshun/ligra).

The adjacency graph format first gives offsets for
verteices, then presents directed edges ordered by their source vertex. The specific format
is as follows:

```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```

### input from SNAP

To transfer any SNAP dataset into a usable format, run the following scripts (here we use DBLP as an example)
```
python3 edgelist_format.py ../../large_graphs/com-dblp.ungraph.txt ../../large_graphs/dblp.txt ../../large_graphs/com-dblp.all.cmty.txt ../../large_graphs/dblp.all
./snap_converter -s -i ../../large_graphs/dblp.txt -o ../../large_graphs/dblp.adj
```
