# motif-based community detection

## Sequential Python implementation

Check the notebook in the sequential_python folder.

## Parallel C++ implementation

The framework is built based on GBBS[https://github.com/ParAlg/gbbs]. The config setup follows GBBS and uses Homegrown scheduler. Compile with bazel[https://bazel.build/install]:
```
cd parallel
bazel build //benchmarks/community/wt:Thresholding_main
```
To run
```
bazel run //benchmarks/community/wt:Thresholding_main -- -rounds 1 -community absolute/path/to/groundtruth/community -method METHOD -threshold DELTA -s absolute/path/to/graph/adj
```
For example, the following command runs TW on DBLP with threshold -15
```
bazel-3.7.0 run //benchmarks/community/wt:Thresholding_main -- -rounds 1 -community path_to_root/large_graphs/dblp.all -method WT -threshold -15.0 -s path_to_root/large_graphs/dblp.adj
```
Available motif-based similarity methods are: TW, TECTONIC, JACCARD and K3.

### input from SNAP

To transfer any SNAP dataset into a usable format, run the following scripts (here we use DBLP as an example)
```
python3 edgelist_format.py ../../large_graphs/com-dblp.ungraph.txt ../../large_graphs/dblp.txt ../../large_graphs/com-dblp.all.cmty.txt ../../large_graphs/dblp.all
./snap_converter -s -i ../../large_graphs/dblp.txt -o ../../large_graphs/dblp.adj
```
