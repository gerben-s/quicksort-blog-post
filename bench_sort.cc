// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>

#include <benchmark/benchmark.h>

#include "hybrid_qsort.h"
#include "third_party/lomuto/lomuto.h"

constexpr int FLAGS_number = 100000;

// Reference impl of the canonical Lomuto partitioning scheme
template <typename T>
size_t LomutoPartition(T* arr, size_t n, T* scratch, size_t scratch_size) {
  auto pivot = arr[n - 1];
  size_t i = 0;
  for (size_t j = 0; j < n - 1; j++) {
    if (arr[j] < pivot) {
      std::swap(arr[i], arr[j]);
      i++;
    }
  }
  std::swap(arr[i], arr[n - 1]);
  return i;
}

template <typename T>
size_t SearchLeaf(T* arr, size_t idx, size_t n) {
  auto left = 2 * idx;
  while (left < n) {
    auto right = left + 1;
    auto l = arr[left];
    auto r = arr[right];
    bool is_smaller = l < r;
    auto max = is_smaller ? r : l;
    arr[idx] = max;
    idx = is_smaller ? right : left;
    left = 2 * idx;
  }
  if (left == n) {
    arr[idx] = arr[left];
    idx = left;
  }
  return idx;
}

template <typename T>
void SiftDown(T* arr, size_t idx, size_t n, T elem) {
  auto left = 2 * idx;
  while (left < n) {
    auto right = left + 1;
    auto l = arr[left];
    auto r = arr[right];
    bool is_smaller = l < r;
    auto max = is_smaller ? r : l;
    auto next_idx = is_smaller ? right : left;
    is_smaller = max < elem;
    arr[idx] = is_smaller ? elem : max;
    elem = is_smaller ? max : elem;
    idx = next_idx;
    left = 2 * next_idx;
  }
  if (left == n) {
    auto l = arr[left];
    bool is_smaller = l < elem;
    arr[idx] = is_smaller ? elem : l;
    elem = is_smaller ? l : elem;
    idx = left;    
  }
  arr[idx] = elem;
}

template <typename T>
void Insert(T* arr, T cur, size_t root, size_t idx) {
  auto parent = idx / 2;
  if (parent >= root) {
    auto tmp = arr[parent];
    if (tmp < cur) {
      arr[idx] = tmp;
      return Insert(arr, cur, root, parent);
    }
  }
  arr[idx] = cur;
}

template <typename T>
void Heapify(T* arr, size_t n) {
  for (int i = n / 2; i >= 1; i--) SiftDown(arr, i, n, arr[i]);
}

template <typename T>
bool IsHeap(T* arr, size_t n) {
  for (size_t i = 1; i <= n / 2; i++) {
    if (arr[i] < arr[2 * i]) return false;
    if (2 * i + 1 <= n && arr[i] < arr[2 * i + 1]) return false;
  }
  return true;
}

template <typename T>
void HeapSort(T* arr, size_t n) {
  arr -= 1;
  Heapify(arr, n);
  for (size_t i = n; i > 1; i--) {
    auto tmp = arr[i];
    arr[i] = arr[1];
#if 0
    SiftDown(arr, 1, i - 1, tmp);
#elif 1
    if (i >= 33) {
      auto idx1 = SearchLeaf(arr, 1, 31);
      if (arr[idx1 >> 1] < tmp) {
        Insert(arr, tmp, 1, idx1);
        continue;
      }
      auto tmp2 = arr[i - 1];
      arr[i - 1] = arr[1];
      auto idx2 = SearchLeaf(arr, 1, 15);
      if (idx2 == idx1 >> 1) {
        // Colliding path
        SiftDown(arr, idx1, i - 2, tmp);
        SiftDown(arr, idx2, i - 2, tmp2);
      } else {
        
      }
    } else {
      SiftDown(arr, 1, i - 1, tmp);
    }
#else
    auto leaf = SearchLeaf(arr, 1, i - 1);
    Insert(arr, tmp, 1, leaf);
#endif
  }
}

template <typename RandomIt, typename Compare>
void InsertionSort(RandomIt first, RandomIt last, Compare comp) {
  auto n = last - first;
  for (decltype(n) i = 1; i < n; i++) {
    auto x = first[i];
    if (comp(x, first[0])) {
      // Move the whole array
      for (auto j = i; j > 0; j--) first[j] = first[j - 1];
      first[0] = x;
    } else {
      // Move part of the array.
      auto j = i;
      for (; comp(x, first[j - 1]); j--) first[j] = first[j - 1];
      first[j] = x;
    }
  }
}

template <typename T, size_t N>
struct Pipe {
  Pipe(const T* arr) {
    for (int i = 0; i < N; i++) pipe_[i] = arr[i];
    exp_gerbens::BubbleSort(pipe_, N);
  }
  T Push(T val) {
    auto res = val;
    for (int i = N - 1; i >= 0; i--) {
      auto tmp = pipe_[i];
      bool smaller = val < tmp;
      pipe_[i] = smaller ? tmp : res;
      res = smaller ? val : tmp;
    }
    return res;
  }
  T pipe_[N];
};

// Generalization of bubble sort. On each iteration the maximum element will
// bubble up the array to it's right place, requiring n iterations of the
// bubbling. If we bubble up the two largest elements we need n / 2 iterations.
// Luckily bubbling up two elements is almost as quick as long as the elements
// can be kept in register.
template <size_t N, typename T>
void PipedBubbleSort(T* arr, size_t n) {
  size_t i;
  for (i = n; i > N - 1; i -= N) {
    Pipe<T, N> pipe(arr);
    for (size_t j = N; j < i; j++) {
      arr[j - N] = pipe.Push(arr[j]);
    }
    for (size_t j = 0; j < N; j++) arr[i - N + j] = pipe.pipe_[j];
  }
  static_assert(N <= 3, "Must patch up remaining elements only supported for N == 3 at the moment");
  if (N == 3 && i == 2) {
    if (arr[1] < arr[0]) std::swap(arr[0], arr[1]);
  }
}

template <typename T>
void MergeFast(const T* left, size_t n, T* out) {
  size_t leftp = 0;
  size_t rightp = 0;
  size_t num_iters = n / 2;
  auto right = left + num_iters;
  for (size_t k = 0; k < num_iters; k++) {
      // If left run head exists and is <= existing right run head.
      auto left1 = left[k - leftp];
      auto left2 = right[leftp];
      bool val = left2 < left1;
      out[k] = val ? left2 : left1;
      leftp += val;
      auto right1 = right[-rightp - 1];
      auto right2 = left[n - 1 - k + rightp];
      val = right2 < right1;
      out[n - 1 - k] = val ? right1 : right2;
      rightp += val;
  }
  if (n & 1) {
    if (num_iters - leftp == num_iters - 1 - rightp) out[num_iters] = left[num_iters - leftp];
    else out[num_iters] = left[n - 1 - num_iters + rightp];
  }
}

void BottomUpMerge(int* arr, size_t n, int* scratch) {
  if (n < 2) return;
  for (size_t i = 0; i < n; i += 2) if (arr[i + 1] < arr[i]) std::swap(arr[i], arr[i + 1]);
  if (n < 4) return;
  for (size_t i = 0; i < n; i += 4) {
    auto a = arr[i];
    auto b = arr[i + 1];
    auto c = arr[i + 2];
    auto d = arr[i + 3];
    arr[i] = c < a ? c : a;
    c = c < a ? a : c;
    arr[i + 3] = d < b ? b : d;
    b = d < b ? d : b;
    arr[i + 1] = b < c ? b : c;
    arr[i + 2] = b < c ? c : b;
  }
  for (size_t w = 8; w <= n; w *= 2) {
    for (size_t i = 0; i < n; i += w) {
      MergeFast(arr + i, w, scratch + i);
    }
    std::swap(arr, scratch);
  }
}

template <typename T>
void Merge(const T* left, size_t n, T* out) {
  auto middle = n / 2;
  size_t smaller_cnt = 0;

  for (size_t k = 0;; k++) {
    auto i = k - smaller_cnt;
    auto j = middle + smaller_cnt;
    if (i >= middle) {
      for (; k < n; k++) out[k] = left[j++];
      return;
    }
    if (j >= n) {
      for (; k < n; k++) out[k] = left[i++];
      return;
    }
    auto x = left[i];
    auto y = left[j];
    bool is_smaller = y < x;
    out[k] = is_smaller ? y : x;
    smaller_cnt += is_smaller;
  }
}


template <bool kDoubleMerge, typename T>
void MergeSort(T* arr, size_t n, T* scratch) {
  if (n <= exp_gerbens::kSmallSortThreshold) return exp_gerbens::BubbleSort2(arr, arr + n, std::less<>{});
  MergeSort<kDoubleMerge>(arr, n / 2, scratch);
  MergeSort<kDoubleMerge>(arr + n / 2, n - n / 2, scratch + n / 2);
  memcpy(scratch, arr, n * sizeof(T));
  if (kDoubleMerge) {
    MergeFast(scratch, n, arr);
  } else {
    Merge(scratch, n, arr);
  }
}


// Benchmarks

template <typename RandomIt>
void std_heap_sort(RandomIt first, RandomIt last) {
  std::make_heap(first, last);
  std::sort_heap(first, last);
}

template <typename RandomIt>
void HeapSort(RandomIt first, RandomIt last) {
  HeapSort(&*first, last - first);
}

template<void (*qsort)(int*, int*)>
void BM_Sort(benchmark::State& state) {
  static std::vector<int> buf(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, 100000000);
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    for (auto& x : buf) x = dist(rnd);
    state.ResumeTiming();
    qsort(buf.data(), buf.data() + FLAGS_number);
  }
}

BENCHMARK_TEMPLATE(BM_Sort, std::sort);
BENCHMARK_TEMPLATE(BM_Sort, std::stable_sort);
BENCHMARK_TEMPLATE(BM_Sort, std_heap_sort);
BENCHMARK_TEMPLATE(BM_Sort, andrei::sort);
BENCHMARK_TEMPLATE(BM_Sort, exp_gerbens::QuickSort);
BENCHMARK_TEMPLATE(BM_Sort, HeapSort);

template<void (*qsort)(int*, int*)>
void BM_SortDuplicates(benchmark::State& state) {
  static std::vector<int> buf(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, 5);
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    for (auto& x : buf) x = dist(rnd);
    state.ResumeTiming();
    qsort(buf.data(), buf.data() + FLAGS_number);
  }
}

BENCHMARK_TEMPLATE(BM_SortDuplicates, std::sort);
BENCHMARK_TEMPLATE(BM_SortDuplicates, std::stable_sort);
BENCHMARK_TEMPLATE(BM_SortDuplicates, std_heap_sort);
BENCHMARK_TEMPLATE(BM_SortDuplicates, andrei::sort);
BENCHMARK_TEMPLATE(BM_SortDuplicates, exp_gerbens::QuickSort);
BENCHMARK_TEMPLATE(BM_SortDuplicates, HeapSort);

template <void (*msort)(int*, size_t, int*)>
void BM_MergeSort(benchmark::State& state) {
  std::vector<int> buf(FLAGS_number);
  std::vector<int> scratch(FLAGS_number);
  std::vector<int> buf3(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(100000000);
  for (auto& x : buf3) x = dist(rnd);
  int power2 = 1;
  while (power2 * 2 <= FLAGS_number) power2 *= 2;
  while (state.KeepRunningBatch(power2)) {
    state.PauseTiming();
    buf = buf3;
    state.ResumeTiming();
    msort(buf.data(), power2, scratch.data());
  }
}
BENCHMARK_TEMPLATE(BM_MergeSort, MergeSort<false>);
BENCHMARK_TEMPLATE(BM_MergeSort, MergeSort<true>);
BENCHMARK_TEMPLATE(BM_MergeSort, BottomUpMerge);

template<void (*small_sort)(int*, int*, std::less<>)>
void BM_SmallSort(benchmark::State& state) {
  int n = state.range(0);
  std::vector<int> buf(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(100000000);
  auto num_batches = FLAGS_number / n * n;
  while (state.KeepRunningBatch(num_batches)) {
    state.PauseTiming();
    for (auto& x : buf) x = dist(rnd);
    state.ResumeTiming();
    for (int i = 0; i + n <= FLAGS_number; i += n) {
      small_sort(buf.data() + i, buf.data() + i + n, std::less<>{});
    }
  }
}
BENCHMARK_TEMPLATE(BM_SmallSort, exp_gerbens::BubbleSort)->Range(2, 32)->RangeMultiplier(2);
BENCHMARK_TEMPLATE(BM_SmallSort, exp_gerbens::BubbleSort2)->Range(2, 32)->RangeMultiplier(2);
BENCHMARK_TEMPLATE(BM_SmallSort, InsertionSort)->Range(2, 32)->RangeMultiplier(2);

void BM_Partition(benchmark::State& state) {
  std::vector<int> buf(FLAGS_number);
  std::vector<int> buf2(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(100000000);
  for (auto& x : buf2) x = dist(rnd);
  constexpr size_t kScratchSize = 512;
  int scratch[kScratchSize];
  auto pivot = buf2[FLAGS_number / 2];
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    buf = buf2;
    state.ResumeTiming();
    auto p = exp_gerbens::HoareLomutoHybridPartition<kScratchSize>(pivot, buf.begin(), buf.end(), scratch, std::less<>{});
    benchmark::DoNotOptimize(p);
  }
}

template <size_t N>
struct PointlessPointerIndirection {
  PointlessPointerIndirection() { me_ = this; }
  PointlessPointerIndirection(const PointlessPointerIndirection&) { me_ = this; }
  PointlessPointerIndirection& operator=(const PointlessPointerIndirection&) {}
  void* Get() {
    auto me = this;
    for (size_t i = 0; i < N; i++) me = me->me_;
    return me;
  }
  PointlessPointerIndirection* me_;
};

template <size_t N>
struct Pointless {
  Pointless() = default;
  Pointless(PointlessPointerIndirection<N>* p) : ptr_(p) {}

  PointlessPointerIndirection<N> *ptr_ = nullptr;
};

template <size_t N>
inline bool operator<(const Pointless<N>& a, const Pointless<N>& b) {
  return a.ptr_->Get() < b.ptr_->Get();
}
template <size_t N>
inline bool operator<=(const Pointless<N>& a, const Pointless<N>& b) {
  return a.ptr_->Get() <= b.ptr_->Get();
}
template <size_t N>
inline bool operator>(const Pointless<N>& a, const Pointless<N>& b) {
  return a.ptr_->Get() > b.ptr_->Get();
}
template <size_t N>
inline bool operator>=(const Pointless<N>& a, const Pointless<N>& b) {
  return a.ptr_->Get() >= b.ptr_->Get();
}


template <size_t N, void (*qsort)(Pointless<N>*, Pointless<N>*)>
void BM_IndirectionSort(benchmark::State& state) {
  std::vector<PointlessPointerIndirection<N>> pointers(FLAGS_number);
  std::vector<Pointless<N>> buf(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number - 1);
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    for (auto& x : buf) x = pointers.data() + dist(rnd);
    state.ResumeTiming();
    qsort(buf.data(), buf.data() + buf.size());
  }
}

BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, std::sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, std::stable_sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, std_heap_sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, andrei::sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, exp_gerbens::QuickSort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, HeapSort);

template <size_t N>
void BM_IndirectionMergeSort(benchmark::State& state) {
  std::vector<PointlessPointerIndirection<N>> pointers(FLAGS_number);
  std::vector<Pointless<N>> buf(FLAGS_number);
  std::vector<Pointless<N>> buf2(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number - 1);
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    for (auto& x : buf) x = pointers.data() + dist(rnd);
    state.ResumeTiming();
    MergeSort<true>(buf.data(), buf.size(), buf2.data());
  }
}

BENCHMARK_TEMPLATE(BM_IndirectionMergeSort, 1);

BENCHMARK_TEMPLATE(BM_IndirectionSort, 0, std::sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 0, std::stable_sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 0, std_heap_sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 0, andrei::sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 0, exp_gerbens::QuickSort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 0, HeapSort);
BENCHMARK_TEMPLATE(BM_IndirectionMergeSort, 0);

template <size_t N>
void BM_IndirectPartition(benchmark::State& state) {
  std::vector<PointlessPointerIndirection<N>> pointers(FLAGS_number);
  std::vector<Pointless<N>> buf(FLAGS_number);
  std::vector<Pointless<N>> buf2(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number - 1);
  for (auto& x : buf2) x = pointers.data() + dist(rnd);
  constexpr size_t kScratchSize = 512;
  Pointless<N> buffer[kScratchSize];
  auto pivot = buf2[FLAGS_number / 2];
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    buf = buf2;
    state.ResumeTiming();
    auto p = exp_gerbens::HoareLomutoHybridPartition<kScratchSize>(pivot, buf.begin(), buf.end(), buffer, std::less<>{});
    benchmark::DoNotOptimize(p);
  }
}
BENCHMARK_TEMPLATE(BM_IndirectPartition, 0);
BENCHMARK_TEMPLATE(BM_IndirectPartition, 1);
BENCHMARK_TEMPLATE(BM_IndirectPartition, 2);
BENCHMARK_TEMPLATE(BM_IndirectPartition, 3);

template <size_t N>
void BM_IndirectMerge(benchmark::State& state) {
  std::vector<PointlessPointerIndirection<N>> pointers(FLAGS_number);
  std::vector<Pointless<N>> buf(FLAGS_number);
  std::vector<Pointless<N>> buf2(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number - 1);
  for (auto& x : buf) x = pointers.data() + dist(rnd);
  while (state.KeepRunningBatch(FLAGS_number)) {
    MergeFast(buf.data(), FLAGS_number, buf2.data());
  }
}
BENCHMARK_TEMPLATE(BM_IndirectMerge, 0);
BENCHMARK_TEMPLATE(BM_IndirectMerge, 1);
BENCHMARK_TEMPLATE(BM_IndirectMerge, 2);
BENCHMARK_TEMPLATE(BM_IndirectMerge, 3);

int main(int argc, char* argv[]) {
  benchmark::Initialize(&argc, argv);

  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number * 4);
  std::vector<int> buf(FLAGS_number);
  for (auto& x : buf) x = dist(rnd);
  auto buf2 = buf;
  auto buf3 = buf;
  auto buf4 = buf;
  // for (auto x : buf) std::cout << x << " "; std::cout << "\n";
  exp_gerbens::QuickSort(buf.begin(), buf.end());
  // for (auto x : buf) std::cout << x << " "; std::cout << "\n";
  std::sort(buf2.begin(), buf2.end());
  // for (auto x : buf2) std::cout << x << " "; std::cout << "\n";
  assert(buf == buf2);
  MergeSort<true>(buf3.data(), buf3.size(), buf2.data());
  assert(buf == buf3);
  HeapSort(buf4.data(), buf4.size());
  assert(buf == buf4);

  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
