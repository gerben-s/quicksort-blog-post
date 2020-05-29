Hoare’s rebuttal and bubble sort’s comeback
Gerben Stavenga and Josh Haberman
Recently Andrei Alexandrescu published an interesting post about optimizing QuickSort using the Lomuto partition scheme. The essence of that post is that for many situations the performance of QuickSort is completely dominated by branch mispredicts and that a big speed up can be achieved by writing branchless code. This has been observed by many and various branchless sorting routines have been proposed. Andrei observed that from the two well known QuickSort partitioning schemes Lomuto is easily implemented branchless, and this indeed performs much better for sorting small primitives. Purely coincidentally we were experimenting with similar ideas but took them in a slightly different direction. We feel that still there is some analysis that complements Andrei’s post and some interesting results are presented.
Basic QuickSort fundamentals
Quicksort refers to a class of algorithms for sorting an array that all share the same outline
void QuickSort(T* left, T* right) {
  if (right - left > kCutOff) {
    auto pivot = ChoosePivotElement(left, right);  // Important but not focus here
    auto p = Partition(pivot, left, right);  // The main work loop
    QuickSort(left, p);
    QuickSort(p, right);  // Tail call, ideally the largest sub-interval
  } else {
    SortSmallArray(left, right);
  }
}
	

Countless variations exist varying in the choice of kCutOff, choice of the sorting algorithm for the small arrays and choice of pivot element. These are important for performance but the main work QuickSort is performing is done in the Partition function. There are two canonical schemes for implementing Partition. The original Hoare scheme and the Lomuto scheme. The Hoare partition scheme works by swapping elements that violate the partition property from the front of the array with elements from the back of the array, processing the array from the outside inwards converging on the partition point somewhere in the middle.
T* HoarePartition(T pivot, T* left, T* right) {
  while (left < right) {
     left = ScanForward(pivot, left, right);
     if (left == right) break;
     right = ScanBackward(pivot, left, right);
     if (left == right) break;
     swap(*left, *right);
  }
  return left;
}
	

In contrast the Lomuto scheme processes the array from front to back, maintaining a properly partitioned array for the elements processed so far at each step.
T* LomutoPartition(T pivot, T* left, T* right) {
  T* p = left;
  for (auto it = left; it < right; it++) {
    if (*it < pivot) {
      std::swap(*it, *p);   // Could be a self-swap
      p++;
    }
  }
  return p;
}
	A remarkably simple loop, with an additional property that it’s stable with respect to the ordering of the elements smaller than the pivot. It’s easy to see that the above loop can be implemented without conditional branches. The conditional swapping can be implemented by conditional moves and the conditional pointer increase could be implemented as a unconditional p += (*it < pivot). Andrei’s blog shows the performance gain of this simple branchless loop over production quality implementations in various standard libraries.
Performance analysis of branchless Lomuto partitioning
Here I want to take the optimizations further and dive deeper into the performance analysis of sorting algorithms. When conditional branching is removed the performance of code tends to become much more stable and easier to understand, data dependencies become key. We are going to analyse the code in terms of a self explanatory pseudo-assembly code, where named variables should be thought of as CPU registers and loads/stores to memory are made explicit. In this notation the basic loop above becomes
loop:
val = Load(it)              // 1
prev_val = Load(p)          // 2
is_smaller = val < pivot    // 3
valnew = cmov(is_smaller, prev_val, val)   // 4
prev_val = cmov(is_smaller, val, prev_val) // 5
Store(it, valnew)           // 6
Store(p, prev_val)          // 7
p += is_smaller             // 8
it++                        // 9
if (it < right) goto loop
	How would it run on a parallel out-of order CPU? The initial model is that the CPU is arbitrary parallel, ie. it can execute as many instructions parallel as long as the instructions are independent. What is important to understand is the loop carried dependency chains. They determine the minimum latency a loop could possibly run. In the above code you see that the dependency between iterations of the loop are carried by it and p. Only line 8 and 9 participate on the loop carried chain and both are single cycle instructions. So we determine that the loop could potentially run at a throughput of 1 cycle per iteration.


  

However things are more complicated in the above loop. There are also loads and stores to memory. If you store something to a certain memory address and later load from that address you must read the value of the previous store. That means loads have dependencies on stores, or at least if the address overlaps. Here it and p are dynamic values and for sure they can overlap, depicted by the dashed lines in the diagram above. So let’s add the fact that there is a dependency between the loads at line 1 and 2 on the stores of line 6 and 7.


  

This completely changes the game, now there is a long loop carried data dependency, lines 6 and 7 depend on lines 4 and 5, which both depend on line 3 , which depend on the loads at lines 1 and 2, which potentially depend on the stores at lines 6 and 7 of the previous iteration. If we count the cycles, we get 5 cycles for the loads (loads itself can be done in parallel), 1 cycle for the comparison, 1 cycle for the conditional move and 1 cycle for the store, hence this loop will run ~8 cycles. A far cry from the 1 cycle iterations our cursory discussion indicated.


Although it’s not possible to reorder stores and loads in general it’s essential for performance of a CPU to do so. Let’s take a simple memcpy loop
loop:
val = Load(it)
Store(it + delta, val)
it++
if (it < end) goto loop
	If the load cannot be reordered with the previous store, this is a 5 + 1 = 6 cycle latency loop. However in memcpy it’s guaranteed that loads and stores never overlap. If the CPU would instead reorder, the above loop would execute with a throughput of one iteration per cycle. It’s execution would look like, ignoring the instructions needed for control flow.
val_0 = Load(it_0); it_1 = it_0 + 1;    // Cycle 1
val_1 = Load(it_1); it_2 = it_1 + 1;    // Cycle 2
val_2 = Load(it_2); it_3 = it_2 + 1;    // Cycle 3
val_3 = Load(it_3); it_4 = it_3 + 1;    // Cycle 4
// The value of the load at cycle 1 becomes available.
// From this point all instructions of the loop are executed each cycle.
val_4 = Load(it_4); it_5 = it_4 + 1; Store(val_0);    // Cycle 5
val_5 = Load(it_5); it_6 = it_5 + 1; Store(val_1);    // Cycle 6
	In practice most of the stores preceding a load in the instruction stream are in fact to different addresses and it is possible to reorder loads in front of stores. Therefore CPU’s do reorder loads in front of stores, which is called speculative loading, however the effects of the load are only committed when it’s verified no store has invalidated the speculative load. If a preceding store, in effect, invalidates the load, execution is rolled back to the load and the CPU starts over. One can imagine that this is very costly and very akin to branch mispredicts. While a lot of stores and loads are to different addresses, there are also plenty of stores and loads to the same address, think about register spills for example. Therefore the CPU uses a prediction model based on instruction address to determine if a load has a dependency on previous stores. In general the CPU is pretty conservative, the cost of a wrong reordering is very high. In the code above the CPU will encounter loads from the same address as recent stores and will be hesitant to do the necessary re-ordering. 
Revisiting the Lomuto partition scheme
Looking closer it’s the load from p that is problematic. The load from it is in fact always from a different address than the previous stores. Furthermore the load at p is also responsible for a lot of extra work in the loop. It is necessary as otherwise the store at p will corrupt values in the array. The values that are overwritten are values previously encountered in the scan and are those elements that are larger than the pivot. If instead we would save these values in a temporary buffer there is no need to swap.
// Distributes the elements [left, right) into begin of left and from end of
// scratch buffer 
T* DistributeForward(T pivot, T* left, T* right, T* scratch) {
  auto scratch_end = scratch + kScratchSize - 1;
  ptrdiff_t offset = 0;
  for (auto it = left; it < right; it++) {
    auto val = *it;
    bool is_larger = val >= pivot;
    auto dst = is_larger ? scratch_end : it;
    dst[offset] = val;
    offset[a][b] -= is_larger;
  }
  return right + offset;
}


T* ModifiedLomutoPartition(T pivot, T* left, T* right, T* scratch) {
  auto p = DistributeForward(pivot, left, right, scratch);
  // To complete the partition we need to copy the elements in scratch
  // to the end of the array.
  auto size = right - p;
  memcpy(p, scratch + kScratchSize - size, size * sizeof(T));
  return p;
}
	This is a much simpler loop, only one load and one store per iteration. More importantly the load will never clash with a previous store. This loop runs much faster than the original loop, it’s not 1 cycle per iteration but 2.5 cycles on my machine. This is indicative that it’s saturating the ILP of the CPU. Unfortunately the above code is not in-place anymore, it requires O(n) additional memory for the scratch buffer.
  

The elegant hybrid
If instead we use a smallish fixed size temporary buffer, we can still use the above code except we need to abort when the fixed buffer is full. What do we do then? It returned the partition point p where it left the loop. At this point [left, p) have the correct elements smaller than pivot at the front of the array. The scratch buffer is full with elements larger or equal to pivot and [p, p + kScratchSize) contains information we don’t need anymore. The idea is that we can do the same algorithm but backwards, we can use [p, p + kScratchSize) as a temporary buffer. Notice how the DistributeForward fills the scratch buffer from back to front, the backwards version would fill the scratch from front to back. So performing DistributeBackwards using the interval  [p, p + kScratchSize) as scratch will neatly pack all smaller elements encountered to the correct place. This continues until the scratch space is full, but now a new scratch space at the end of the array opened up. Wait, this looks like Hoare’s algorithm but hybridized with the lomuto inspired distribute function. 
T* ModifiedHoarePartition(T pivot, T* left, T* right, T* scratch) {
  auto pleft = DistributeForward(left, right, pivot, scratch);
  if (right - pleft <= kScratchSize) {
    auto size = right - pleft;
    std::memcpy(pleft, scratch + kScratchSize - size, size * sizeof(T));
    return pleft;
  }
  left = pleft + kScratchSize;
  T* res;
  while (true) {
    right = DistributeBackward(left, right, pivot, left - kScratchSize)
        - kScratchSize;
    if (right <= left) {
      res = right;
      break;
    }
    left = DistributeForward(left, right, pivot, right) + kScratchSize;
    if (right <= left) {
      res = left - kScratchSize;
      break;
    }
  }
  std::memcpy(res, scratch, kScratchSize * sizeof(T));
  return res;
}
	What we are ending up with is an in-place algorithm that’s almost branchfree. The precise number of iterations in the modified Lomuto partitioning scheme depend on the outcome of the comparisons and will be difficult to predict exactly right. However the loop is guaranteed to iterate at least kScratchSize times. This basically amortises the cost of branch mispredictions over many elements making them irrelevant for performance. I consider this to be a truly elegant design.
Fallback for sorting short arrays
The next step is the fallback for short arrays where the overhead of recursion starts dominating. In the literature insertion sort is most often recommended together with countless micro optimizations applied. I found that at this point partitioning was so fast that QuickSort beat insertion sort all the way down to just a few elements. The problem is that insertion sort has unpredictable branches, basically 1 miss per insert. The solution is Bubble sort. Bubble sort has a very predictable access pattern and the swap can be implemented branchless. A little more optimizing you discover you don’t need a swap. One can keep the maximum in register and store the minimum.
void BubbleSort(T* arr, size_t n) {
  for (size_t i = n; i > 1; i--) {
    auto max = arr[0];
    for (size_t j = 1; j < i; j++) {
      auto y = arr[j];
      arr[j - 1] = (max <= y ? max : y);
      max = (max <= y ? y : max);
    }
    arr[i - 1] = max;
  }
}
	In the above inner-loop max is the variable that participates on the longest loop carried data chain, with a compare and conditional move on it. This makes the above loop execute in 2 cycles per iteration. However we can do better. If instead of bubbling the max, we’re bubbling up the two largest elements we only need to iterate the bubbling stage n/2 times, instead of n times. It turns out that using a clever implementation one can bubble two elements in the same 2 cycles. In fact it’s possible to generalize this, one can bubble up the top N in constant cycles per iteration in the model where the CPU has arbitrary ILP.
void BubbleSort2(T* arr, size_t n) {
  for (size_t i = n; i > 1; i -= 2) {
    auto x = arr[0];
    auto y = arr[1];
    if (y < x) std::swap(x, y);
    for (size_t j = 2; j < i; j++) {
      auto z = arr[j];
      bool is_smaller = y <= z;
      auto w = is_smaller ? y : z;
      y = is_smaller ? z : y;
      is_smaller = x <= z
      arr[j - 2] = (is_smaller ? x : z);
      x = is_smaller ? w : x;
    }
    arr[i - 2] = x;
    arr[i - 1] = y;
  }
}
	The benchmarks verify that indeed this makes a difference.
BM_SmallSort<exp_gerbens::BubbleSort>/2                1.63           1.63   429410000  
BM_SmallSort<exp_gerbens::BubbleSort>/8                3.15           3.15   222330000  
BM_SmallSort<exp_gerbens::BubbleSort>/32              10.4           10.4     67322112  
BM_SmallSort<exp_gerbens::BubbleSort2>/2               1.36           1.36   514920000  
BM_SmallSort<exp_gerbens::BubbleSort2>/8               2.18           2.18   320980000  
BM_SmallSort<exp_gerbens::BubbleSort2>/32              6.71           6.71   100009728 
Bringing it all together
What are the results? The results are nothing short of spectacular. The following is a simple benchmark on sorting 100000 int’s. The time is normalized by the number of elements, so it’s the amount of time spent per element. I’m using clang + libc++ here, gcc is dramatically worse in emitting branch free code.


CPU: Intel Skylake Xeon with HyperThreading (36 cores) dL1:32KB dL2:1024KB dL3:24MB
Benchmark                         Time(ns)        CPU(ns)     Iterations
------------------------------------------------------------------------
BM_Sort<std_sort>                       51.6           51.6     10000000  
BM_Sort<std_stable_sort>                65.6           65.6     10000000  
BM_Sort<lib_qsort>                      90.4           90.5      7800000  
BM_Sort<andrei_sort>                    32.6           32.6     21500000  
BM_Sort<exp_gerbens::QuickSort>         16.4           16.4     43200000 


We’re talking about a 2x win over Andrei’s implementation as copied from his github page.


We’ve seen how crucial it is to understand data dependencies in order to optimize code. Especially hidden memory dependencies between load and stores can greatly influence performance of work loops. Understanding the data dependency graph of code is often where the real performance gains lie, yet very little attention is given to it in the blogosphere. I’ve read many articles about the impact of branch mispredictions, importance of data locality and caches, but much less about data dependencies. I bet that a question like “why are linked lists slow?” is answered by many in terms of locality, caches or unpredictable random memory access. At least I’ve heard those reasons often, even Stroustrup says as much. Those reasons can play a part, but it’s not the main reason. Fundamentally iterating a linked list has a load-to-use on the critical path, making it 5 times slower than iterating a flat array. Furthermore accessing flat arrays allow loop unrolling which can further improve ILP. 
Why is QuickSort fast compared to merge and heap sort?
This brings us to answer why QuickSort is fast compared to other sorts with good or even better theoretical complexity. It’s all about data dependencies. The quick sort partition loop above demonstrates a distinctive feature. The element it will process next does not depend on the outcome of the comparisons in previous iterations. Compare this to merge sort. In merge sort the two head elements are compared, however the next elements that need to be compared depend on the outcome of this comparison. It’s trivial to implement merge sort branch free. It will look like
val1 = Load(left + k - right_idx)  // 1
val2 = Load(right + right_idx)
is_smaller = val2 < val1   // 2
tmp = cmov(is_smaller, val2, val1)
Store(out + k, tmp);
k++;
right_idx += is_smaller  // 3
	This is about 8 cycles per iteration to update right_idx, we have a load and non-trivial indexing at line 1 (6 cycles), and line 2 and 3 both being 1 cycle.  Similar analysis holds for heap sort, restoring the heap property requires comparing the two children and recurse on the subtree of the biggest child. 
left = Load(2 * idx)
right = Load(2 * idx + 1)
is_smaller = left < right
tmp = cmov(is_smaller, right, left)
Store(idx, tmp);
idx = 2 * idx + is_smaller
	Again this is 8 cycles on the loop carried data chain. This goes to the heart of the matter of why QuickSort is fast even though theoretically it has inferior behavior compared to heap and merge sort. By construction heap and merge sort divide the data very evenly in the implied tree structure, the heap structure is a binary tree of minimum depth as is the recursive subdivision that merge sort performs. This means that the number of comparisons they do is which tightly hugs the information theoretical lower bound of . Instead QuickSort bets on obtaining a reasonably balanced tree with high probability. The bits of information extracted from a comparison with the pivot depends on its rank in the array, only pivot’s that are close to the median will result in obtaining 1 bit of information per comparison. Solving the recurrence equation of QuickSort with a uniform pivot choice, gives as the number of comparisons on average. Hence QuickSort does a factor more comparisons or equivalently more iterations in the inner work loop. However the enormous speed difference in the basic work loop more than compensates for this information theoretical factor. Also, when partitioning big arrays spending a little amount of work improving the choice of pivot by a median of three brings this factor down to ~1.2, further improvements can get this factor rather close to 1. The main point here is that this factor is dominated by the differences in throughput of the work loop.






We can significantly speed up merge sort by this one simple trick. Due to the large dependency chain merge sort loop runs at a very low IPC. This basically means we can add more instructions for free. In particular merge sort has a backwards equivalent. We can merge forward and backward in a single loop while keeping the latency of the loop body the same. It also eliminates an awkward exit condition as now you can unconditionally iterate n/2 times. This reduces the number of iterations roughly by 2x. Another trick is preloading values for the next iteration,
next_val1 = Load(left + k - right_idx + 1)
next_val2 = Load(right + right_idx + 1)
is_smaller = val2 < val1
tmp = cmov(is_smaller, val2, val1)
val1 = cmov(is_smaller, val1, next_val1)
val2 = cmov(is_smaller, next_val2, val2)
Store(out + k, tmp);
k++;
right_idx += is_smaller
	You see that in a single iteration the increase of right_idx does not depend on a load as val1 and val2 are already available at the start of the iteration. Over two iterations one can see that right_idx depends on itself. This chain is ~8 cycles long but spread over 2 iterations which gives a throughput of ~4 cycles per iteration. Combining these two tricks could lead to merge sort on par with the simple partition loop of QuickSort. However it’s just a mitigation. If instead of sorting a simple primitive value we would sort pointers to a struct. The comparison operator would have an extra load which immediately adds to the critical path. QuickSort is immune to this, even rather costly comparisons do not influence the ability to make progress. Of course if the comparison itself suffers from frequent branch misses, then that will limit the ability to overlap different stages of the iteration in the CPU pipeline. 




[a]Could offset be increasing instead? I'm guessing you chose this way to save a sub instruction for dst[offset], but that's off the critical path, right? Having offset negative makes things more confusing. For example I think your return statement is wrong, shouldn't it be return + offset?
[b]Yes it turns out that the critical path is just 1 cycle and this loop is fully dominated by instruction count in it's throughput