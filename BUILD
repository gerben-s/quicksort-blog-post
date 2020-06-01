cc_binary(
    name = "bench_sort",
    testonly = 1,
    srcs = [
        "bench_sort.cc",
        "hybrid_qsort.h",
        "third_party/lomuto/lomuto.h",
    ],
    deps = [
        "@com_github_google_benchmark//:benchmark_main",
    ],
    copts = ["-std=c++17"],
)
