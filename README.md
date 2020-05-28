# Branchfree QuickSort performance analysis

This directory contains code and benchmarks accompanying the blog post at

The `main()` function validates the correctness of all the algorithms
before running benchmarks.

The benchmark requires
[Bazel](https://github.com/bazelbuild/bazel/releases) to build.
At build time Bazel will automatically download and build the
benchmarking framework.

From this directory, run:

```
$ bazel build -c opt :benchmark
$ ../bazel-bin/quicksort-blog-post/benchmark
```

