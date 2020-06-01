#include <cassert>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <random>
#include <algorithm>

using std::swap;

using TYPE = long;
static const size_t SORT_THRESHOLD = 16;

/**
Partitions the range [first, last) around a pivot chosen as the minimum of
first[0] and last[-1]. Uses the Hoare partition algorithm.
Returns: a pointer to the new position of the pivot.
*/
long* hoare_partition(long* first, long* last) {
    assert(first <= last);
    if (last - first < 2)
        return first; // nothing interesting to do
    --last;
    if (*first > *last)
        swap(*first, *last);
    auto pivot_pos = first;
    auto pivot = *pivot_pos;
    for (;;) {
        ++first;
        auto f = *first;
        while (f < pivot)
            f = *++first;
        auto l = *last;
        while (pivot < l)
            l = *--last;
        if (first >= last)
            break;
        *first = l;
        *last = f;
        --last;
    }
    --first;
    swap(*first, *pivot_pos);
    return first;
}

/**
Partitions the range [first, last) around a pivot chosen as the minimum of
first[0] and last[-1]. Uses the Lomuto partition algorithm.
Returns: a pointer to the new position of the pivot.
*/
long* lomuto_partition(long* first, long* last) {
    assert(first <= last);
    if (last - first < 2)
        return first; // nothing interesting to do
    // Choose pivot.
    --last;
    if (*first > *last)
        swap(*first, *last);
    auto pivot_pos = first;
    auto pivot = *first;
    // Prelude: position first (aka the write head) right on the first larger element.
    do {
        ++first;
    } while (*first < pivot);
    assert(first <= last);
    // Main loop.
    for (auto read = first + 1; read < last; ++read) {
        auto x = *read;
        if (x < pivot) {
            *read = *first;
            *first = x;
            ++first;
        }
    }
    // Move the pivot to its final slot.
    assert(*first >= pivot);
    --first;
    *pivot_pos = *first;
    *first = pivot;
    return first;
}

/**
Partitions the range [first, last) around a pivot chosen as the minimum of
first[0] and last[-1]. Uses the Lomuto partition algorithm, branch-free.
Returns: a pointer to the new position of the pivot.
*/
long* lomuto_partition_branchfree(long* first, long* last) {
    assert(first <= last);
    if (last - first < 2)
        return first; // nothing interesting to do
    // Choose pivot.
    --last;
    if (*first > *last)
        swap(*first, *last);
    auto pivot_pos = first;
    auto pivot = *first;
    // Prelude.
    do {
        ++first;
        assert(first <= last);
    } while (*first < pivot);
    // Main loop.
    for (auto read = first + 1; read < last; ++read) {
        auto x = *read;
        auto less = -int(x < pivot);
        auto delta = less & (read - first);
        first[delta] = *first;
        read[-delta] = x;
        first -= less;
    }
    // Move the pivot to its final slot.
    assert(*first >= pivot);
    --first;
    *pivot_pos = *first;
    *first = pivot;
    return first;
}

/**
Inserts *last into the range to its left, assumed to be sorted and have at
least one element smaller than *last.
*/
template<typename It>
void unguarded_linear_insert(It last) {
    auto val = *last;
    --last;
    auto x = *last;
    if (val >= x)
        return;
    for (;;) {
        last[1] = x;
        --last;
        x = *last;
        if (val >= x)
            break;
    }
    last[1] = val;
}

/**
Sorts [first, last) using insertion sort.
*/
template<typename It>
void insertion_sort(It first, It last) {
    assert(first <= last);
    for (auto i = first + 1; i < last; ++i) {
        auto val = *i;
	    if (val < *first) {
            size_t n = i - first - 1;
            do {
                first[n + 1] = first[n];
            }
            while (n--);
	        *first = val;
	    }
	    else
	        unguarded_linear_insert(i);
	}
}

/**
Sorts [first, last) using quicksort and insertion sort for short subarrays.
*/
template <class It>
void sort(It first, It last) {
    while (last - first > size_t(SORT_THRESHOLD)) {
#if defined(LOMUTO_BRANCHY)
	    auto cut = lomuto_partition(first, last);
#elif defined(LOMUTO)
	    auto cut = lomuto_partition_branchfree(first, last);
#else
	    auto cut = hoare_partition(first, last);
#endif
        assert(cut >= first);
        assert(cut < last);
	    sort(cut + 1, last);
	    last = cut;
	}
    insertion_sort(first, last);
}

/**
Returns the difference between two timespecs, in milliseconds.
*/
double timediff(const timespec& start, const timespec& end) {
    return (end.tv_nsec - start.tv_nsec) / double(1e6)
        + (end.tv_sec - start.tv_sec) * double(1e3);
}

/**
Verifies that an array has v[i]==i for all i.
*/
void checkData(std::vector<TYPE>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == i) continue;
        fprintf(stderr, "Array has been corrupted at position %zu.\n", i);
        abort();
    }
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    const size_t length = atol(argv[1]);
    if (length == 0) return 2;
    const size_t repeats = 50000000 / length;
    std::vector<double> times(repeats);
    timespec start, end;
    double totalTime = 0, minTime = std::numeric_limits<double>::max();
    std::mt19937 g(1942);
    std::vector<TYPE> v(length);
    
    for (size_t i = 0; i < v.size(); ++i) v[i] = i;

    for (size_t i = 0; i < repeats; ++i) {
        std::shuffle(v.begin(), v.end(), g);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
#if defined(STDSORT)
        std::sort(&v.front(), 1 + &v.back());
#else
        sort(&v.front(), 1 + &v.back());
#endif
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        // Make sure everything is correct
        checkData(v);
        // Bookkeeping
        auto t = timediff(start, end);
        times[i] = t;
        if (t < minTime) {
            minTime = t;
        }
        totalTime += t;
    }

    std::sort(times.begin(), times.end());
    printf("min_milliseconds=%.4f\n", minTime);
    printf("median_milliseconds=%.4f\n", times[times.size() / 2]);
}
