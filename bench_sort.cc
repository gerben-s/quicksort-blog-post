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
size_t SiftDown(T* arr, size_t idx, size_t n) {
  auto left = 2 * idx;
  auto right = 2 * idx + 1;
  if (right <= n) {
    auto l = arr[left];
    auto r = arr[right];
    bool is_smaller = l < r;
    auto max = is_smaller ? r : l;
    auto next_idx = is_smaller ? right : left;
    arr[idx] = max;
    return SiftDown(arr, next_idx, n);
  } else if (left == n) {
    arr[idx] = arr[left];
    return left;
  }
  return idx;
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
  for (int i = n / 2; i >= 1; i--) {
    auto tmp = arr[i];
    auto j = SiftDown(arr, i, n);
    Insert(arr, tmp, i, j);
  }
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
    auto tmp = arr[1];
    auto j = SiftDown(arr, 1, i);
    Insert(arr, arr[i], 1, j);
    arr[i] = tmp;
  }
}

// OddEven has no inter iteration dependencies but a lot of stores.
// #stores if N/2 * (N + N) = N^2
template <typename T>
void OddEvenSort(T* arr, size_t n) {
  for (size_t i = 0; i < (n + 1) / 2; i++) {
    for (size_t j = 0; j < n - 1; j += 2) exp_gerbens::BranchlessSwap(arr, j, j + 1);
    for (size_t j = 1; j < n - 1; j += 2) exp_gerbens::BranchlessSwap(arr, j, j + 1);
  }
}

template <typename T, size_t N>
struct Pipe {
  Pipe(const T* arr) {
    static_assert(N > 0 && N <= 3, "");
    switch (N) {
      case 1:
        pipe_[0] = arr[0];
        break;
      case 2: {
        auto p0 = arr[0];
        auto p1 = arr[1];
        if (p1 < p0) std::swap(p0, p1);
        pipe_[0] = p0;
        pipe_[1] = p1;
        break;
      }
      case 3: {
        auto p0 = arr[0];
        auto p1 = arr[1];
        auto p2 = arr[2];
        if (p1 < p0) std::swap(p0, p1);
        if (p2 < p1) std::swap(p1, p2);
        if (p1 < p0) std::swap(p0, p1);
        pipe_[0] = p0;
        pipe_[1] = p1;
        pipe_[2] = p2;
        break;
      }
    }
  }
  T Push(T val) {
    bool smaller = val < pipe_[0];
    auto res = smaller ? val : pipe_[0];
    pipe_[0] = smaller ? pipe_[0] : val;
    for (int i = 1; i < N; i++) {
      smaller = val < pipe_[i];
      pipe_[i - 1] = smaller ? pipe_[i - 1] : pipe_[i];
      pipe_[i] = smaller ? pipe_[i] : val;
//      asm("" : "+r"(pipe_[i]));
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
  if (i == 2) {
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
  auto iMiddle = n / 2;
  size_t  i = 0, j = iMiddle;

  // While there are elements in the left or right runs...
  for (size_t k = 0; k < n; k++) {
    // If left run head exists and is <= existing right run head.
    if (i < iMiddle && (j >= n || left[i] <= left[j])) {
        out[k] = left[i];
        i = i + 1;
    } else {
        out[k] = left[j];
        j = j + 1;
    }
  }
}

template <typename T>
void MergeSort(T* arr, size_t n, T* scratch) {
  if (n <= exp_gerbens::kSmallSortThreshold) return exp_gerbens::BubbleSort2(arr, n);
  MergeSort(arr, n / 2, scratch);
  MergeSort(arr + n / 2, n - n / 2, scratch + n / 2);
  memcpy(scratch, arr, n * sizeof(T));
  constexpr bool kFastMerge = true;
  if (kFastMerge) {
    MergeFast(scratch, n, arr);
  } else {
    Merge(scratch, n, arr);
  }
}


// Benchmarks

template <typename T>
void std_sort(T* x, size_t n) { std::sort(x, x + n); }
template <typename T>
void std_heap_sort(T* x, size_t n) {
  std::make_heap(x, x + n);
  std::sort_heap(x, x + n);
}
template <typename T>
void std_stable_sort(T* x, size_t n) { std::stable_sort(x, x + n); }

void lib_qsort(int* x, size_t n) {
  std::qsort(x, n, 4, [](const void* a, const void* b) { return *static_cast<const int*>(a) - *static_cast<const int*>(b); });
}

template<void (*qsort)(int*, size_t)>
void BM_Sort(benchmark::State& state) {
  static std::vector<int> buf(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(100000000);
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    for (auto& x : buf) x = dist(rnd);
    state.ResumeTiming();
    qsort(buf.data(), buf.size());
  }
}

BENCHMARK_TEMPLATE(BM_Sort, std_sort);
BENCHMARK_TEMPLATE(BM_Sort, std_stable_sort);
BENCHMARK_TEMPLATE(BM_Sort, std_heap_sort);
BENCHMARK_TEMPLATE(BM_Sort, lib_qsort);
BENCHMARK_TEMPLATE(BM_Sort, exp_gerbens::QuickSort);
BENCHMARK_TEMPLATE(BM_Sort, HeapSort);

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
BENCHMARK_TEMPLATE(BM_MergeSort, MergeSort);
BENCHMARK_TEMPLATE(BM_MergeSort, BottomUpMerge);

template<void (*small_sort)(int*, size_t)>
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
      small_sort(buf.data() + i, n);
    }
  }
}
BENCHMARK_TEMPLATE(BM_SmallSort, OddEvenSort)->Range(2, 32);
BENCHMARK_TEMPLATE(BM_SmallSort, exp_gerbens::BubbleSort)->Range(2, 32)->RangeMultiplier(2);
BENCHMARK_TEMPLATE(BM_SmallSort, exp_gerbens::BubbleSort2)->Range(2, 32)->RangeMultiplier(2);

void BM_Partition(benchmark::State& state) {
  std::vector<int> buf(FLAGS_number);
  std::vector<int> buf2(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(100000000);
  for (auto& x : buf2) x = dist(rnd);
  constexpr size_t kScratchSize = 512;
  int scratch[kScratchSize];
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    buf = buf2;
    state.ResumeTiming();
    auto pivot = exp_gerbens::Partition<kScratchSize>(buf.data(), FLAGS_number, scratch);
    benchmark::DoNotOptimize(pivot);
  }
}
BENCHMARK(BM_Partition);
void BM_PartitionInto(benchmark::State& state) {
  std::vector<int> buf(FLAGS_number);
  std::vector<int> buf2(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(100000000);
  for (auto& x : buf) x = dist(rnd);
  while (state.KeepRunningBatch(FLAGS_number)) {
    exp_gerbens::PartitionInto(buf.data(), FLAGS_number, buf2.data());
  }
}
BENCHMARK(BM_PartitionInto);

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


template <size_t N, void (*qsort)(Pointless<N>*, size_t)>
void BM_IndirectionSort(benchmark::State& state) {
  std::vector<PointlessPointerIndirection<N>> pointers(FLAGS_number);
  std::vector<Pointless<N>> buf(FLAGS_number);
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number - 1);
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    for (auto& x : buf) x = pointers.data() + dist(rnd);
    state.ResumeTiming();
    qsort(buf.data(), buf.size());
  }
}

BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, std_sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, std_stable_sort);
BENCHMARK_TEMPLATE(BM_IndirectionSort, 1, std_heap_sort);
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
    MergeSort(buf.data(), buf.size(), buf2.data());
  }
}

BENCHMARK_TEMPLATE(BM_IndirectionMergeSort, 1);

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
  while (state.KeepRunningBatch(FLAGS_number)) {
    state.PauseTiming();
    buf = buf2;
    state.ResumeTiming();
    auto p = exp_gerbens::Partition<kScratchSize>(buf.data(), FLAGS_number, buffer);
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
  std::mt19937 rnd;
  std::uniform_int_distribution<int> dist(0, FLAGS_number * 4);
  std::vector<int> buf(FLAGS_number);
  for (auto& x : buf) x = dist(rnd);
  auto buf2 = buf;
  auto buf3 = buf;
  auto buf4 = buf;
  // for (auto x : buf) std::cout << x << " "; std::cout << "\n";
  exp_gerbens::QuickSort(buf.data(), buf.size());
  // for (auto x : buf) std::cout << x << " "; std::cout << "\n";
  std::sort(buf2.begin(), buf2.end());
  // for (auto x : buf2) std::cout << x << " "; std::cout << "\n";
  assert(buf == buf2);
  MergeSort(buf3.data(), buf3.size(), buf2.data());
  assert(buf == buf3);
  HeapSort(buf4.data(), buf4.size());
  assert(buf == buf4);

  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
