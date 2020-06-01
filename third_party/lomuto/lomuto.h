#include <cassert>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <random>
#include <algorithm>

namespace andrei {

using std::swap;

using TYPE = long;
static const size_t SORT_THRESHOLD = 16;

/**
Partitions the range [first, last) around a pivot chosen as the minimum of
first[0] and last[-1]. Uses the Hoare partition algorithm.
Returns: a pointer to the new position of the pivot.
*/
template <typename T>
T* hoare_partition(T* first, T* last) {
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
template <typename T>
T* lomuto_partition(T* first, T* last) {
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
template <typename T>
T* lomuto_partition_branchfree(T* first, T* last) {
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
	    auto cut = lomuto_partition_branchfree(first, last);
        assert(cut >= first);
        assert(cut < last);
	    sort(cut + 1, last);
	    last = cut;
	}
    insertion_sort(first, last);
}

}  // namespace andrei
