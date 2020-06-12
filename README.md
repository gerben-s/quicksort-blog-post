# Branchfree QuickSort performance analysis

This is not an officially supported Google product

This directory contains code and benchmarks accompanying the blog post about QuickSort performance.

The `main()` function validates the correctness of the main algorithms before running benchmarks.

The benchmark requires [Bazel](https://github.com/bazelbuild/bazel/releases) to build.
At build time Bazel will automatically download and build the benchmarking framework.

From this directory, run:

```
$ CC=clang bazel build -c opt :bench_sort
$ ../bazel-bin/quicksort-blog-post/bench_sort
```

The versions of gcc I've run did not lower conditionals into branchfree code.

# Overview of benchmarks

Random int's 
```
Benchmark                                          Time           CPU Iterations
---------------------------------------------------------------------------------
BM_Sort<std::sort>                                79 ns         78 ns    9300000
BM_Sort<std::stable_sort>                         90 ns         90 ns    7400000
BM_Sort<std_heap_sort>                           130 ns        130 ns    4600000
BM_Sort<andrei::sort>                             52 ns         52 ns   15100000
BM_Sort<exp_gerbens::QuickSort>                   30 ns         30 ns   24700000
BM_Sort<pdqsort>                                  42 ns         42 ns   16800000
BM_Sort<HeapSort>                                 51 ns         51 ns   14500000
```
Random pointers sort on address (0 levels of indirection)
```
Benchmark                                              Time           CPU Iterations
---------------------------------------------------------------------------------
BM_IndirectionSort<0, std::sort>                      77 ns         77 ns    9200000
BM_IndirectionSort<0, std::stable_sort>               92 ns         91 ns    7600000
BM_IndirectionSort<0, std_heap_sort>                 124 ns        124 ns    5800000
BM_IndirectionSort<0, andrei::sort>                   56 ns         56 ns   10000000
BM_IndirectionSort<0, exp_gerbens::QuickSort>         32 ns         32 ns   18300000
BM_IndirectionSort<0, pdqsort_branchless>             40 ns         40 ns   17600000
BM_IndirectionSort<0, HeapSort>                       60 ns         60 ns   11900000
```
Random pointers sort on value pointed to (1 levels of indirection)
```
Benchmark                                              Time           CPU Iterations
-------------------------------------------------------------------------------------
BM_IndirectionSort<1, std::sort>                      97 ns         97 ns    7400000
BM_IndirectionSort<1, std::stable_sort>              133 ns        133 ns    5100000
BM_IndirectionSort<1, std_heap_sort>                 180 ns        180 ns    4100000
BM_IndirectionSort<1, andrei::sort>                   67 ns         67 ns   11600000
BM_IndirectionSort<1, exp_gerbens::QuickSort>         42 ns         42 ns   16300000
BM_IndirectionSort<1, pdqsort_branchless>             54 ns         54 ns   10000000
BM_IndirectionSort<1, HeapSort>                      131 ns        131 ns    6000000
```
