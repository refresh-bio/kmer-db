#pragma once
// *****************************************************************************
// This header file is based on Pattern-defeating quicksort implementation by Orson Peters.
// It extends the original code by adding parallel variants of the algorithms.
// 
// -----------------------------------------------------------------------------
// The original licence by Orson Peters:
// -----------------------------------------------------------------------------
// pdqsort.h - Pattern - defeating quicksort.
// 
// Copyright(c) 2021 Orson Peters
// 
// This software is provided 'as-is', without any express or implied warranty.In no event will the
// authors be held liable for any damages arising from the use of this software.
// 
// Permission is granted to anyone to use this software for any purpose, including commercial
// applications, and to alter it and redistribute it freely, subject to the following restrictions :
// 
// 1. The origin of this software must not be misrepresented; you must not claim that you wrote the
// original software.If you use this software in a product, an acknowledgment in the product
// documentation would be appreciated but is not required.
// 
// 2. Altered source versions must be plainly marked as such, and must not be misrepresented as
// being the original software.
// 
// 3. This notice may not be removed or altered from any source distribution.
// -----------------------------------------------------------------------------
// *****************************************************************************

#include <future>
#include <vector>
#include <mutex>
#include <queue>
#include <tuple>
#include <functional>
#include <map>
#include <list>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <type_traits>
#include <cinttypes>

#ifdef _WIN32
#include <emmintrin.h>
#endif

#include "../../active_thread_pool/lib/active_thread_pool.h"

#define REFRESH_USE_THREAD_POOLS

//#define REFRESH_ATOMIC_QUEUE_VECTOR
#define REFRESH_ATOMIC_QUEUE_MAP
//#define REFRESH_ATOMIC_QUEUE_LIST

namespace refresh
{
    // *****************************************************************************
    //
    // *****************************************************************************
    namespace sort
    {
        namespace details
        {
            template<typename T>
            class atomic_queue
            {
                std::atomic_flag a_lock;
    #if defined(REFRESH_ATOMIC_QUEUE_VECTOR)
                std::vector<T> data_queue;
    #elif defined(REFRESH_ATOMIC_QUEUE_MAP)
                std::multimap<ptrdiff_t, T> data_queue;
    #elif defined(REFRESH_ATOMIC_QUEUE_LIST)
                std::queue<T, std::list<T>> data_queue;
    #else
                std::queue<T> data_queue;
    #endif

                std::atomic_flag is_completed;
                size_t n_pause;

                inline void pause(size_t n = 1)
                {
                    for (size_t i = 0; i < n; ++i)
    #if defined(_WIN32)
                        _mm_pause();
    #elif defined(__aarch64__)
                        std::this_thread::yield();
    #else
                        __builtin_ia32_pause();
    #endif
                }

                inline void push_lock_acquire()
                {
                    while (a_lock.test_and_set(std::memory_order_acquire))
                        while (a_lock.test(std::memory_order_relaxed))
                            ;// pause();
                }

                inline void pop_lock_acquire()
                {
                    while (a_lock.test_and_set(std::memory_order_acquire))
                        while (a_lock.test(std::memory_order_relaxed))
                            pause(n_pause);
                }

                inline void lock_clear()
                {
                    a_lock.clear(std::memory_order_release);
                }

                void adjust_n_pause(size_t n_threads)
                {
                    n_pause = 2;
                    while (n_threads >>= 1)
                        ++n_pause;
                }

            public:
                atomic_queue(size_t n_threads = 1)
                {
    #if defined(REFRESH_ATOMIC_QUEUE_VECTOR)
                    data_queue.reserve(16 + n_threads);
    #endif
                    adjust_n_pause(n_threads);
                }

                void push(const T &new_value, const size_t priority = 0)
                {
                    push_lock_acquire();
    #if defined(REFRESH_ATOMIC_QUEUE_VECTOR)
                    data_queue.push_back(move(new_value));
    #elif defined(REFRESH_ATOMIC_QUEUE_MAP)
                    data_queue.emplace(priority, move(new_value));
    #else
                    data_queue.emplace(move(new_value));
    #endif

                    lock_clear();
                }

                bool wait_and_pop(T& value)
                {
                    while (true)
                    {
                        if (is_completed.test())
                            return false;

                        pop_lock_acquire();

                        if (!data_queue.empty())
                            break;

                        lock_clear();
                        pause(n_pause);
                    }

    #if defined(REFRESH_ATOMIC_QUEUE_VECTOR)
                    value = move(data_queue.back());
                    data_queue.pop_back();
    #elif defined(REFRESH_ATOMIC_QUEUE_MAP)
                    value = move(data_queue.rbegin()->second);
                    data_queue.erase(--data_queue.end());
    #else
                    value = move(data_queue.front());
                    data_queue.pop();
    #endif

                    lock_clear();

                    return true;
                }

                void restart(size_t n_threads = 1)
                {
                    push_lock_acquire();

                    is_completed.clear();
                    adjust_n_pause(n_threads);

                    lock_clear();
                }

                void set_completed()
                {
                    is_completed.test_and_set();
                }
            };

            template<class T> struct is_default_compare : std::false_type { };
            template<class T> struct is_default_compare<std::less<T>> : std::true_type { };
            template<class T> struct is_default_compare<std::greater<T>> : std::true_type { };

            // Returns floor(log2(n)), assumes n > 0.
            template<class T>
            static inline int log2(T n) {
                int log = 0;
                while (n >>= 1) ++log;
                return log;
            }
        }

        // *****************************************************************************
        //
        // *****************************************************************************
        template<typename Iter, typename Compare, typename ThreadPool, bool Branchless>
        class pdq_sorter
        {
            Compare comp;
            ThreadPool& thread_pool;

            using job_t = std::tuple<Iter, Iter, int, bool>;

            details::atomic_queue<job_t> jobs;

#ifdef REFRESH_USE_THREAD_POOLS
            ThreadPool::pool_state_t pool_state;
#else
            std::vector<std::future<void>> workers;
#endif

            std::atomic<int> active_threads{};
            std::atomic_flag threads_ready{};

            std::atomic<size_t> no_started_workers = 0;
            size_t max_started_workers = 0;

            void join_workers()
            {
#ifdef REFRESH_USE_THREAD_POOLS
                pool_state.heavy_wait();
#else
                for (size_t i = 0; i < std::min<size_t>(max_started_workers, no_started_workers); ++i)
                    workers[i].wait();

                workers.clear();
#endif
            }

            void add_worker(size_t id)
            {
#ifdef REFRESH_USE_THREAD_POOLS
                thread_pool.launch([&]
                    {
                        job_t job;

                        bool is_active = thread_pool.is_active();

                        while (true)
                        {
                            if (!jobs.wait_and_pop(job))
                                break;

                            pdqsort_loop_par(std::get<0>(job), std::get<1>(job), comp, std::get<2>(job), std::get<3>(job), is_active);
                            if (!decrease_active_threads())
                                break;
                        }
                    },
                    & pool_state);
#else
                workers[id] = std::async(std::launch::async, [&]
                    {
                        job_t job;

                        while (true)
                        {
                            if (!jobs.wait_and_pop(job))
                                break;

                            pdqsort_loop_par(std::get<0>(job), std::get<1>(job), comp, std::get<2>(job), std::get<3>(job), false);
                            if (!decrease_active_threads())
                                break;
                        }
                    });
#endif
            }

            void increase_active_threads()
            {
                ++active_threads;
            }

            bool decrease_active_threads()
            {
                auto x = active_threads.fetch_sub(1);

                if (x == 1)
                {
                    jobs.set_completed();

//                    threads_ready.test_and_set();
//                    threads_ready.notify_one();

                    return false;
                }

                return true;
            }

            // Partitions below this size are sorted using insertion sort.
            static constexpr size_t insertion_sort_threshold = 24;        // 24

            // Partitions above this size use Tukey's ninther to select the pivot.
            static constexpr size_t ninther_threshold = 128;              // 128

            // When we detect an already sorted partition, attempt an insertion sort that allows this
            // amount of element moves before giving up.
            static constexpr size_t partial_insertion_sort_limit = 8;

            // Must be multiple of 8 due to loop unrolling, and < 256 to fit in unsigned char.
            static constexpr size_t block_size = 64;

            // Cacheline size, assumes power of two.
            static constexpr size_t cacheline_size = 64;

            // Min job size for parallel processing
//            static constexpr size_t min_par_job_size = 1024;
            static constexpr size_t min_par_job_size[2] = { 1024, 768 };

            // Sorts [begin, end) using insertion sort with the given comparison function.
            static inline void insertion_sort(Iter begin, Iter end, Compare comp) {
                typedef typename std::iterator_traits<Iter>::value_type T;
                if (begin == end) return;

                for (Iter cur = begin + 1; cur != end; ++cur) {
                    Iter sift = cur;
                    Iter sift_1 = cur - 1;

                    // Compare first so we can avoid 2 moves for an element already positioned correctly.
                    if (comp(*sift, *sift_1)) {
                        T tmp = std::move(*sift);

                        do { *sift-- = std::move(*sift_1); } while (sift != begin && comp(tmp, *--sift_1));

                        *sift = std::move(tmp);
                    }
                }
            }

            // Sorts [begin, end) using insertion sort with the given comparison function. Assumes
            // *(begin - 1) is an element smaller than or equal to any element in [begin, end).
            static inline void unguarded_insertion_sort(Iter begin, Iter end, Compare comp) {
                typedef typename std::iterator_traits<Iter>::value_type T;
                if (begin == end) return;

                for (Iter cur = begin + 1; cur != end; ++cur) {
                    Iter sift = cur;
                    Iter sift_1 = cur - 1;

                    // Compare first so we can avoid 2 moves for an element already positioned correctly.
                    if (comp(*sift, *sift_1)) {
                        T tmp = std::move(*sift);

                        do { *sift-- = std::move(*sift_1); } while (comp(tmp, *--sift_1));

                        *sift = std::move(tmp);
                    }
                }
            }

            // Attempts to use insertion sort on [begin, end). Will return false if more than
            // partial_insertion_sort_limit elements were moved, and abort sorting. Otherwise it will
            // successfully sort and return true.
            static inline bool partial_insertion_sort(Iter begin, Iter end, Compare comp) {
                typedef typename std::iterator_traits<Iter>::value_type T;
                if (begin == end) return true;

                std::size_t limit = 0;
                for (Iter cur = begin + 1; cur != end; ++cur) {
                    Iter sift = cur;
                    Iter sift_1 = cur - 1;

                    // Compare first so we can avoid 2 moves for an element already positioned correctly.
                    if (comp(*sift, *sift_1)) {
                        T tmp = std::move(*sift);

                        do { *sift-- = std::move(*sift_1); } while (sift != begin && comp(tmp, *--sift_1));

                        *sift = std::move(tmp);
                        limit += cur - sift;
                    }

                    if (limit > partial_insertion_sort_limit) return false;
                }

                return true;
            }

            static inline void sort2(Iter a, Iter b, Compare comp) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            }

            // Sorts the elements *a, *b and *c using comparison function comp.
            static inline void sort3(Iter a, Iter b, Iter c, Compare comp) {
                sort2(a, b, comp);
                sort2(b, c, comp);
                sort2(a, b, comp);

                /*            auto x0 = std::move(*a);
                            auto x1 = std::move(*b);
                            auto x2 = std::move(*c);

                            int o0 = 0, o1 = 1, o2 = 0;
                            auto r0 = comp(x1, x0);	o0 += r0;	o1 -= r0;
                            auto r1 = comp(x2, x0);	o0 += r1;
                            r0 = comp(x2, x1);	o1 += r0;
                            o2 = 3 - o0 - o1;

                            Iter arr[3] = { a, b, c };

                            *(arr[o0]) = std::move(x0);
                            *(arr[o1]) = std::move(x1);
                            *(arr[o2]) = std::move(x2);*/
            }

            template<class T>
            static inline T* align_cacheline(T* p) {
                std::uintptr_t ip = reinterpret_cast<std::uintptr_t>(p);
                while (ip % 64 != 0)
                    ++ip;
                return reinterpret_cast<T*>(ip);
            }

            static inline void swap_offsets(Iter first, Iter last,
                unsigned char* offsets_l, unsigned char* offsets_r,
                size_t num, bool use_swaps) {
                typedef typename std::iterator_traits<Iter>::value_type T;
                if (use_swaps) {
                    // This case is needed for the descending distribution, where we need
                    // to have proper swapping for pdqsort to remain O(n).
                    for (size_t i = 0; i < num; ++i) {
                        std::iter_swap(first + offsets_l[i], last - offsets_r[i]);
                    }
                }
                else if (num > 0) {
                    Iter l = first + offsets_l[0]; Iter r = last - offsets_r[0];
                    T tmp(std::move(*l)); *l = std::move(*r);
                    for (size_t i = 1; i < num; ++i) {
                        l = first + offsets_l[i]; *r = std::move(*l);
                        r = last - offsets_r[i]; *l = std::move(*r);
                    }
                    *r = std::move(tmp);
                }
            }

            // Partitions [begin, end) around pivot *begin using comparison function comp. Elements equal
            // to the pivot are put in the right-hand partition. Returns the position of the pivot after
            // partitioning and whether the passed sequence already was correctly partitioned. Assumes the
            // pivot is a median of at least 3 elements and that [begin, end) is at least
            // insertion_sort_threshold long. Uses branchless partitioning.
            static inline std::pair<Iter, bool> partition_right_branchless(Iter begin, Iter end, Compare comp) {
                typedef typename std::iterator_traits<Iter>::value_type T;

                // Move pivot into local for speed.
                T pivot(std::move(*begin));
                Iter first = begin;
                Iter last = end;

                // Find the first element greater than or equal than the pivot (the median of 3 guarantees
                // this exists).
                while (comp(*++first, pivot));

                // Find the first element strictly smaller than the pivot. We have to guard this search if
                // there was no element before *first.
                if (first - 1 == begin) while (first < last && !comp(*--last, pivot));
                else                    while (!comp(*--last, pivot));

                // If the first pair of elements that should be swapped to partition are the same element,
                // the passed in sequence already was correctly partitioned.
                bool already_partitioned = first >= last;
                if (!already_partitioned) {
                    std::iter_swap(first, last);
                    ++first;

                    // The following branchless partitioning is derived from "BlockQuicksort: How Branch
                    // Mispredictions don’t affect Quicksort" by Stefan Edelkamp and Armin Weiss, but
                    // heavily micro-optimized.
                    unsigned char offsets_l_storage[block_size + cacheline_size];
                    unsigned char offsets_r_storage[block_size + cacheline_size];
                    unsigned char* offsets_l = align_cacheline(offsets_l_storage);
                    unsigned char* offsets_r = align_cacheline(offsets_r_storage);

                    Iter offsets_l_base = first;
                    Iter offsets_r_base = last;
                    size_t num_l, num_r, start_l, start_r;
                    num_l = num_r = start_l = start_r = 0;

                    while (first < last) {
                        // Fill up offset blocks with elements that are on the wrong side.
                        // First we determine how much elements are considered for each offset block.
                        size_t num_unknown = last - first;
                        size_t left_split = num_l == 0 ? (num_r == 0 ? num_unknown / 2 : num_unknown) : 0;
                        size_t right_split = num_r == 0 ? (num_unknown - left_split) : 0;

                        // Fill the offset blocks.
                        if (left_split >= block_size) [[likely]] {
                            for (size_t i = 0; i < block_size;) {
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                            }
                            }
                        else {
                            for (size_t i = 0; i < left_split;) {
                                offsets_l[num_l] = (unsigned char)i++; num_l += !comp(*first, pivot); ++first;
                            }
                        }

                        if (right_split >= block_size) [[likely]] {
                            for (size_t i = 0; i < block_size;) {
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                            }
                            }
                        else {
                            for (size_t i = 0; i < right_split;) {
                                offsets_r[num_r] = (unsigned char) ++i; num_r += comp(*--last, pivot);
                            }
                        }

                        // Swap elements and update block sizes and first/last boundaries.
                        size_t num = std::min(num_l, num_r);
                        swap_offsets(offsets_l_base, offsets_r_base,
                            offsets_l + start_l, offsets_r + start_r,
                            num, num_l == num_r);
                        num_l -= num; num_r -= num;
                        start_l += num; start_r += num;

                        if (num_l == 0) {
                            start_l = 0;
                            offsets_l_base = first;
                        }

                        if (num_r == 0) {
                            start_r = 0;
                            offsets_r_base = last;
                        }
                    }

                    // We have now fully identified [first, last)'s proper position. Swap the last elements.
                    if (num_l) {
                        offsets_l += start_l;
                        while (num_l--) std::iter_swap(offsets_l_base + offsets_l[num_l], --last);
                        first = last;
                    }
                    if (num_r) {
                        offsets_r += start_r;
                        while (num_r--) std::iter_swap(offsets_r_base - offsets_r[num_r], first), ++first;
                        last = first;
                    }
                }

                // Put the pivot in the right place.
                Iter pivot_pos = first - 1;
                *begin = std::move(*pivot_pos);
                *pivot_pos = std::move(pivot);

                return std::make_pair(pivot_pos, already_partitioned);
            }

            // Partitions [begin, end) around pivot *begin using comparison function comp. Elements equal
            // to the pivot are put in the right-hand partition. Returns the position of the pivot after
            // partitioning and whether the passed sequence already was correctly partitioned. Assumes the
            // pivot is a median of at least 3 elements and that [begin, end) is at least
            // insertion_sort_threshold long.
            static inline std::pair<Iter, bool> partition_right(Iter begin, Iter end, Compare comp) {
                typedef typename std::iterator_traits<Iter>::value_type T;

                // Move pivot into local for speed.
                T pivot(std::move(*begin));

                Iter first = begin;
                Iter last = end;

                // Find the first element greater than or equal than the pivot (the median of 3 guarantees
                // this exists).
                while (comp(*++first, pivot));

                // Find the first element strictly smaller than the pivot. We have to guard this search if
                // there was no element before *first.
                if (first - 1 == begin) while (first < last && !comp(*--last, pivot));
                else                    while (!comp(*--last, pivot));

                // If the first pair of elements that should be swapped to partition are the same element,
                // the passed in sequence already was correctly partitioned.
                bool already_partitioned = first >= last;

                // Keep swapping pairs of elements that are on the wrong side of the pivot. Previously
                // swapped pairs guard the searches, which is why the first iteration is special-cased
                // above.
                while (first < last) {
                    std::iter_swap(first, last);
                    while (comp(*++first, pivot));
                    while (!comp(*--last, pivot));
                }

                // Put the pivot in the right place.
                Iter pivot_pos = first - 1;
                *begin = std::move(*pivot_pos);
                *pivot_pos = std::move(pivot);

                return std::make_pair(pivot_pos, already_partitioned);
            }

            // Similar function to the one above, except elements equal to the pivot are put to the left of
            // the pivot and it doesn't check or return if the passed sequence already was partitioned.
            // Since this is rarely used (the many equal case), and in that case pdqsort already has O(n)
            // performance, no block quicksort is applied here for simplicity.
            static inline Iter partition_left(Iter begin, Iter end, Compare comp) {
                typedef typename std::iterator_traits<Iter>::value_type T;

                T pivot(std::move(*begin));
                Iter first = begin;
                Iter last = end;

                while (comp(pivot, *--last));

                if (last + 1 == end) while (first < last && !comp(pivot, *++first));
                else                 while (!comp(pivot, *++first));

                while (first < last) {
                    std::iter_swap(first, last);
                    while (comp(pivot, *--last));
                    while (!comp(pivot, *++first));
                }

                Iter pivot_pos = last;
                *begin = std::move(*pivot_pos);
                *pivot_pos = std::move(pivot);

                return pivot_pos;
            }

            static inline void pdqsort_loop(Iter begin, Iter end, Compare comp, int bad_allowed, bool leftmost, bool is_active) {
                typedef typename std::iterator_traits<Iter>::difference_type diff_t;

                // Use a while loop for tail recursion elimination.
                while (true) {
                    diff_t size = end - begin;

                    // Insertion sort is faster for small arrays.
                    if (size < (diff_t) insertion_sort_threshold) {
                        if (leftmost) insertion_sort(begin, end, comp);
                        else unguarded_insertion_sort(begin, end, comp);
                        return;
                    }

                    // Choose pivot as median of 3 or pseudomedian of 9.
                    diff_t s2 = size / 2;
                    if (size > (diff_t) ninther_threshold) {
                        sort3(begin, begin + s2, end - 1, comp);
                        sort3(begin + 1, begin + (s2 - 1), end - 2, comp);
                        sort3(begin + 2, begin + (s2 + 1), end - 3, comp);
                        sort3(begin + (s2 - 1), begin + s2, begin + (s2 + 1), comp);
                        std::iter_swap(begin, begin + s2);
                    }
                    else sort3(begin + s2, begin, end - 1, comp);

                    // If *(begin - 1) is the end of the right partition of a previous partition operation
                    // there is no element in [begin, end) that is smaller than *(begin - 1). Then if our
                    // pivot compares equal to *(begin - 1) we change strategy, putting equal elements in
                    // the left partition, greater elements in the right partition. We do not have to
                    // recurse on the left partition, since it's sorted (all equal).
                    if (!leftmost && !comp(*(begin - 1), *begin)) {
                        begin = partition_left(begin, end, comp) + 1;
                        continue;
                    }

                    // Partition and get results.
                    std::pair<Iter, bool> part_result =
                        Branchless ? partition_right_branchless(begin, end, comp)
                        : partition_right(begin, end, comp);
                    Iter pivot_pos = part_result.first;
                    bool already_partitioned = part_result.second;

                    // Check for a highly unbalanced partition.
                    diff_t l_size = pivot_pos - begin;
                    diff_t r_size = end - (pivot_pos + 1);
                    bool highly_unbalanced = l_size < size / 8 || r_size < size / 8;

                    // If we got a highly unbalanced partition we shuffle elements to break many patterns.
                    if (highly_unbalanced) [[unlikely]] {
                        // If we had too many bad partitions, switch to heapsort to guarantee O(n log n).
                        if (--bad_allowed == 0) [[unlikely]] {
                            std::make_heap(begin, end, comp);
                            std::sort_heap(begin, end, comp);
                            return;
                            }

                            if (l_size >= (diff_t)insertion_sort_threshold) {
                                std::iter_swap(begin, begin + l_size / 4);
                                std::iter_swap(pivot_pos - 1, pivot_pos - l_size / 4);

                                if (l_size > (diff_t) ninther_threshold) {
                                    std::iter_swap(begin + 1, begin + (l_size / 4 + 1));
                                    std::iter_swap(begin + 2, begin + (l_size / 4 + 2));
                                    std::iter_swap(pivot_pos - 2, pivot_pos - (l_size / 4 + 1));
                                    std::iter_swap(pivot_pos - 3, pivot_pos - (l_size / 4 + 2));
                                }
                            }

                            if (r_size >= (diff_t) insertion_sort_threshold) {
                                std::iter_swap(pivot_pos + 1, pivot_pos + (1 + r_size / 4));
                                std::iter_swap(end - 1, end - r_size / 4);

                                if (r_size > (diff_t) ninther_threshold) {
                                    std::iter_swap(pivot_pos + 2, pivot_pos + (2 + r_size / 4));
                                    std::iter_swap(pivot_pos + 3, pivot_pos + (3 + r_size / 4));
                                    std::iter_swap(end - 2, end - (1 + r_size / 4));
                                    std::iter_swap(end - 3, end - (2 + r_size / 4));
                                }
                            }
                        }
                    else {
                        // If we were decently balanced and we tried to sort an already partitioned
                        // sequence try to use insertion sort.
                        if (already_partitioned && partial_insertion_sort(begin, pivot_pos, comp)
                            && partial_insertion_sort(pivot_pos + 1, end, comp)) return;
                    }

                    // Sort the left partition first using recursion and do tail recursion elimination for
                    // the right-hand partition.
                    pdqsort_loop(begin, pivot_pos, comp, bad_allowed, leftmost, is_active);
                    begin = pivot_pos + 1;
                    leftmost = false;
                }
            }

            inline void pdqsort_loop_par(Iter begin, Iter end, Compare comp, int bad_allowed, bool leftmost, bool is_active) {
                typedef typename std::iterator_traits<Iter>::difference_type diff_t;

                // Use a while loop for tail recursion elimination.
                while (true) {
                    diff_t size = end - begin;

                    // Insertion sort is faster for small arrays.
                    if (size < (diff_t)insertion_sort_threshold) {
                        if (leftmost) insertion_sort(begin, end, comp);
                        else unguarded_insertion_sort(begin, end, comp);
                        return;
                    }

                    // Choose pivot as median of 3 or pseudomedian of 9.
                    diff_t s2 = size / 2;
                    if (size > (diff_t) ninther_threshold) {
                        sort3(begin, begin + s2, end - 1, comp);
                        sort3(begin + 1, begin + (s2 - 1), end - 2, comp);
                        sort3(begin + 2, begin + (s2 + 1), end - 3, comp);
                        sort3(begin + (s2 - 1), begin + s2, begin + (s2 + 1), comp);
                        std::iter_swap(begin, begin + s2);
                    }
                    else sort3(begin + s2, begin, end - 1, comp);

                    // If *(begin - 1) is the end of the right partition of a previous partition operation
                    // there is no element in [begin, end) that is smaller than *(begin - 1). Then if our
                    // pivot compares equal to *(begin - 1) we change strategy, putting equal elements in
                    // the left partition, greater elements in the right partition. We do not have to
                    // recurse on the left partition, since it's sorted (all equal).
                    if (!leftmost && !comp(*(begin - 1), *begin)) {
                        begin = partition_left(begin, end, comp) + 1;
                        continue;
                    }

                    // Partition and get results.
                    std::pair<Iter, bool> part_result =
                        Branchless ? partition_right_branchless(begin, end, comp)
                        : partition_right(begin, end, comp);
                    Iter pivot_pos = part_result.first;
                    bool already_partitioned = part_result.second;

                    // Check for a highly unbalanced partition.
                    diff_t l_size = pivot_pos - begin;
                    diff_t r_size = end - (pivot_pos + 1);
                    bool highly_unbalanced = l_size < size / 8 || r_size < size / 8;

                    // If we got a highly unbalanced partition we shuffle elements to break many patterns.
                    if (highly_unbalanced) [[unlikely]] {
                        // If we had too many bad partitions, switch to heapsort to guarantee O(n log n).
                        if (--bad_allowed == 0) [[unlikely]] {
                            std::make_heap(begin, end, comp);
                            std::sort_heap(begin, end, comp);
                            return;
                            }

                            if (l_size >= (diff_t)insertion_sort_threshold) {
                                std::iter_swap(begin, begin + l_size / 4);
                                std::iter_swap(pivot_pos - 1, pivot_pos - l_size / 4);

                                if (l_size > (diff_t)ninther_threshold) {
                                    std::iter_swap(begin + 1, begin + (l_size / 4 + 1));
                                    std::iter_swap(begin + 2, begin + (l_size / 4 + 2));
                                    std::iter_swap(pivot_pos - 2, pivot_pos - (l_size / 4 + 1));
                                    std::iter_swap(pivot_pos - 3, pivot_pos - (l_size / 4 + 2));
                                }
                            }

                            if (r_size >= (diff_t)insertion_sort_threshold) {
                                std::iter_swap(pivot_pos + 1, pivot_pos + (1 + r_size / 4));
                                std::iter_swap(end - 1, end - r_size / 4);

                                if (r_size > (diff_t)ninther_threshold) {
                                    std::iter_swap(pivot_pos + 2, pivot_pos + (2 + r_size / 4));
                                    std::iter_swap(pivot_pos + 3, pivot_pos + (3 + r_size / 4));
                                    std::iter_swap(end - 2, end - (1 + r_size / 4));
                                    std::iter_swap(end - 3, end - (2 + r_size / 4));
                                }
                            }
                        }
                    else {
                        // If we were decently balanced and we tried to sort an already partitioned
                        // sequence try to use insertion sort.
                        if (already_partitioned && partial_insertion_sort(begin, pivot_pos, comp)
                            && partial_insertion_sort(pivot_pos + 1, end, comp)) return;
                    }

                    // Sort the left partition first using recursion and do tail recursion elimination for
                    // the right-hand partition.

    //                if (size >= min_par_job_size && jobs.size() <= 4 * max_started_workers)
                    if (size >= (diff_t)min_par_job_size[(int) is_active])
                    {
                        increase_active_threads();
                        size_t priority = pivot_pos - begin;
                        if (priority & 1)
                            priority = ~priority;

                        jobs.push(std::make_tuple(begin, pivot_pos, bad_allowed, leftmost), priority);

                        if (size >= (2 - (int) is_active) * (diff_t)min_par_job_size[(int) is_active])
                        {
                            auto worker_id = no_started_workers.fetch_add(1);
                            if (worker_id < max_started_workers)
                                add_worker(worker_id);
                        }
                    }
                    else
                        pdqsort_loop(begin, pivot_pos, comp, bad_allowed, leftmost, is_active);

                    begin = pivot_pos + 1;
                    leftmost = false;
                }
            }

        public:
            pdq_sorter(const Compare& comp, ThreadPool &thread_pool) : 
                comp(comp),
                thread_pool(thread_pool)
            {}

            ~pdq_sorter()
            {
 //               jobs.set_completed();
            }

            static void local_adjust_n_threads(size_t& n_threads, size_t no_items, bool is_active)
            {
                if (n_threads == 0)
                    n_threads = std::thread::hardware_concurrency();

                if (n_threads * min_par_job_size[(int) is_active] > no_items)
                    n_threads = no_items / min_par_job_size[(int) is_active];

                if (n_threads == 0)
                    n_threads = 1;
            }

            static void sort(Iter begin, Iter end, Compare comp)
            {
                auto bad_allowed = details::log2(end - begin);

//                cerr << "\t\t"s + std::to_string(std::distance(begin, end)) + "                    \r";
//                fflush(stderr);

                pdqsort_loop(begin, end, comp, bad_allowed, true, false);
            }

            void sort(size_t n_threads, Iter begin, Iter end)
            {
//                local_adjust_n_threads(n_threads, std::distance(begin, end));

//                cerr << "\t\t"s + std::to_string(std::distance(begin, end)) + "                    \r";
//                fflush(stderr);

                    
                if (n_threads == 1)
                {
                    sort(begin, end, comp);
                    return;
                }

                active_threads = 1;
                auto bad_allowed = details::log2(end - begin);

                jobs.restart(n_threads - 1);

                //            workers.clear();
#ifdef REFRESH_USE_THREAD_POOLS
                pool_state.reserve(n_threads - 1);
#else
                workers.resize(n_threads - 1);
#endif
                max_started_workers = n_threads - 1;
                no_started_workers = 0;

                bool is_active_tp = thread_pool.is_active();

                pdqsort_loop_par(begin, end, comp, bad_allowed, true, is_active_tp);
                decrease_active_threads();

                job_t job;

                while (true)
                {
                    if (!jobs.wait_and_pop(job))
                        break;

                    pdqsort_loop_par(std::get<0>(job), std::get<1>(job), comp, std::get<2>(job), std::get<3>(job), is_active_tp);
                    decrease_active_threads();
                }

//                threads_ready.wait(false);

                join_workers();
            }
        };

        // *****************************************************************************
        // 
        // *****************************************************************************
        inline size_t pdqsort_adjust_threads(size_t size, size_t max_n_threads)
        {
            size_t proposed_n_threads;

            if (size <= 1024)
                proposed_n_threads = 1;
            else if (size <= 2 * 1024)
                proposed_n_threads = 2;
            else if (size <= 4 * 1024)
                proposed_n_threads = 3;
            else if (size <= 12 * 1024)
                proposed_n_threads = 5;
            else if (size <= 16 * 1024)
                proposed_n_threads = 6;
            else if (size <= 32 * 1024)
                proposed_n_threads = 8;
            else if (size <= 64 * 1024)
                proposed_n_threads = 10;
            else if (size <= 128 * 1024)
                proposed_n_threads = 12;
            else if (size <= 256 * 1024)
                proposed_n_threads = 14;
            else if (size <= 512 * 1024)
                proposed_n_threads = 16;
            else if (size <= 768 * 1024)
                proposed_n_threads = 18;
            else if (size <= 1024 * 1024)
                proposed_n_threads = 20;
            else if (size <= 2 * 1024 * 1024)
                proposed_n_threads = 24;
            else if (size <= 4 * 1024 * 1024)
                proposed_n_threads = 28;
            else if (size <= 8 * 1024 * 1024)
                proposed_n_threads = 32;
            else if (size <= 16 * 1024 * 1024)
                proposed_n_threads = 36;
            else if (size <= 32 * 1024 * 1024)
                proposed_n_threads = 40;
            else if (size <= 64 * 1024 * 1024)
                proposed_n_threads = 48;
            else
                proposed_n_threads = 64;

            return std::min<size_t>(proposed_n_threads, max_n_threads);
        }

        // *****************************************************************************
        template<class Iter, class Compare>
        inline void pdqsort(Iter begin, Iter end, Compare comp)
        {
            if (begin == end)
                return;

            pdq_sorter<Iter, Compare, typename refresh::async_pool,
                details::is_default_compare<typename std::decay<Compare>::type>::value&&
                std::is_arithmetic<typename std::iterator_traits<Iter>::value_type>::value>::sort(begin, end, comp);
        }

        // *****************************************************************************
        template<class Iter>
        inline void pdqsort(Iter begin, Iter end)
        {
            typedef typename std::iterator_traits<Iter>::value_type T;
            pdqsort(begin, end, std::less<T>());
        }

        // *****************************************************************************
        template<class Iter, class Compare>
        inline void pdqsort(size_t n_threads, Iter begin, Iter end, Compare comp)
        {
            refresh::async_pool ap;
            pdsqort_tp(n_threads, begin, end, comp, ap);
        }
            
        // *****************************************************************************
        template<class Iter, class Compare, class ThreadPool>
        inline void pdqsort_tp(size_t n_threads, Iter begin, Iter end, Compare comp, ThreadPool &thread_pool)
        {
            pdq_sorter<Iter, Compare, ThreadPool, details::is_default_compare<typename std::decay<Compare>::type>::value&&
                std::is_arithmetic<typename std::iterator_traits<Iter>::value_type>::value>::local_adjust_n_threads(n_threads, std::distance(begin, end), thread_pool.is_active());

            if (n_threads == 1)
                pdqsort<Iter, Compare>(begin, end, comp);
            else
            {
                if (begin == end)
                    return;

                pdq_sorter<Iter, Compare, ThreadPool,
                    details::is_default_compare<typename std::decay<Compare>::type>::value&&
                    std::is_arithmetic<typename std::iterator_traits<Iter>::value_type>::value>(comp, thread_pool).sort(n_threads, begin, end);
            }
        }

        // *****************************************************************************
        template<class Iter, class ThreadPool>
        inline void pdqsort_tp(size_t n_threads, Iter begin, Iter end, ThreadPool &thread_pool)
        {
            typedef typename std::iterator_traits<Iter>::value_type T;
            pdqsort_tp(n_threads, begin, end, std::less<T>(), thread_pool);
        }

        // *****************************************************************************
        template<class Iter>
        inline void pdqsort(size_t n_threads, Iter begin, Iter end)
        {
            refresh::async_pool ap;
            pdqsort_tp(n_threads, begin, end, ap);
        }

        // *****************************************************************************
        template<class Iter, class Compare>
        inline void pdqsort_branchless(Iter begin, Iter end, Compare comp)
        {
            if (begin == end)
                return;

            pdq_sorter<Iter, Compare, typename refresh::async_pool, true>::sort(begin, end, comp);
        }

        // *****************************************************************************
        template<class Iter>
        inline void pdqsort_branchless(Iter begin, Iter end)
        {
            typedef typename std::iterator_traits<Iter>::value_type T;
            pdqsort_branchless(begin, end, std::less<T>());
        }

        // *****************************************************************************
        template<class Iter, class Compare, class ThreadPool>
        inline void pdqsort_branchless_tp(size_t n_threads, Iter begin, Iter end, Compare comp, ThreadPool &thread_pool)
        {
            pdq_sorter<Iter, Compare, ThreadPool, true>::local_adjust_n_threads(n_threads, std::distance(begin, end), thread_pool.is_active());

            if (n_threads == 1)
                pdqsort_branchless(begin, end, comp);
            else
            {
                if (begin == end)
                    return;

                pdq_sorter<Iter, Compare, ThreadPool, true>(comp, thread_pool).sort(n_threads, begin, end);
            }
        }

        // *****************************************************************************
        template<class Iter, class Compare>
        inline void pdqsort_branchless(size_t n_threads, Iter begin, Iter end, Compare comp)
        {
            refresh::async_pool ap;
            pdqsort_branchless_tp(n_threads, begin, end, comp, ap);
        }

        // *****************************************************************************
        template<class Iter, class ThreadPool>
        inline void pdqsort_branchless_tp(size_t n_threads, Iter begin, Iter end, ThreadPool &thread_pool)
        {
            typedef typename std::iterator_traits<Iter>::value_type T;
            pdqsort_branchless_tp(n_threads, begin, end, std::less<T>(), thread_pool);
        }

        // *****************************************************************************
        template<class Iter>
        inline void pdqsort_branchless(size_t n_threads, Iter begin, Iter end)
        {
            refresh::async_pool ap;
            pdqsort_branchless_tp(n_threads, begin, end, ap);
        }
    }
}

//#undef REFRESH_ATOMIC_QUEUE_VECTOR
#undef REFRESH_ATOMIC_QUEUE_MAP
//#undef REFRESH_ATOMIC_QUEUE_LIST

// EOF
