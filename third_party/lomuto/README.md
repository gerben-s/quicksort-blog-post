This code accompanies the article ["Lomuto's Comeback"](https://dlang.org/blog/2020/05/14/lomutos-comeback/) in The D Blog.

To build the benchmarks, use for gdc:

    # std::sort baseline
    gdc -O3 -frelease -fno-assert -fbounds-check=off -fno-postconditions -fno-preconditions -fversion=STDSORT lomuto.d
    # Use Hoare partition
    gdc -O3 -frelease -fno-assert -fbounds-check=off -fno-postconditions -fno-preconditions lomuto.d
    # Use Lomuto partition, traditional implementation
    gdc -O3 -frelease -fno-assert -fbounds-check=off -fno-postconditions -fno-preconditions -fversion=LOMUTO_BRANCHY lomuto.d
    # Use Lomuto partition, branch-free implementation
    gdc -O3 -frelease -fno-assert -fbounds-check=off -fno-postconditions -fno-preconditions -fversion=LOMUTO lomuto.d

For ldc:

    # std::sort baseline
    ldc 1.17.0 -O3 -boundscheck=off -release -d-version=STDSORT lomuto.d
    # Use Hoare partition
    ldc 1.17.0 -O3 -boundscheck=off -release lomuto.d
    # Use Lomuto partition, traditional implementation
    ldc 1.17.0 -O3 -boundscheck=off -release -d-version=LOMUTO_BRANCHY lomuto.d
    # Use Lomuto partition, branch-free implementation
    ldc 1.17.0 -O3 -boundscheck=off -release -d-version=LOMUTO lomuto.d

For gcc:

    # std::sort baseline
    g++ -std=c++17 -O3 -DNDEBUG -DSTDSORT lomuto.cpp
    # Use Hoare partition
    g++ -std=c++17 -O3 -DNDEBUG lomuto.cpp
    # Use Lomuto partition, traditional implementation
    g++ -std=c++17 -O3 -DNDEBUG -DLOMUTO_BRANCHY lomuto.cpp
    # Use Lomuto partition, branch-free implementation
    g++ -std=c++17 -O3 -DNDEBUG -DLOMUTO lomuto.cpp

For clang:

    # std::sort baseline
    clang++ -std=c++17 -O3 -DNDEBUG -DSTDSORT lomuto.cpp
    # Use Hoare partition
    clang++ -std=c++17 -O3 -DNDEBUG lomuto.cpp
    # Use Lomuto partition, traditional implementation
    clang++ -std=c++17 -O3 -DNDEBUG -DLOMUTO_BRANCHY lomuto.cpp
    # Use Lomuto partition, branch-free implementation
    clang++ -std=c++17 -O3 -DNDEBUG -DLOMUTO lomuto.cpp

Run the produced binary with one numeric argument, the number of items to sort. Beware: all compilers write to `./a.out` except for ldc, which writes to `./lomuto`.
