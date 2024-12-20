#ifndef PQUEUE_HPP
#define PQUEUE_HPP

#include "dsGlobal.hpp"
#include <ext/pb_ds/priority_queue.hpp>

MOLSQ_NAMESPACE_HPP_START

template<typename T, typename Cmp = less<T> >
using PQueue = __gnu_pbds::priority_queue<T, Cmp>;

template<typename T, typename Cmp = less<T> >
using BinaryHeap = __gnu_pbds::priority_queue<T, Cmp, __gnu_pbds::binary_heap_tag>;

template<typename T, typename Cmp = less<T> >
using BinomialHeap = __gnu_pbds::priority_queue<T, Cmp, __gnu_pbds::binomial_heap_tag>;

template<typename T, typename Cmp = less<T> >
using RCBinomialHeap = __gnu_pbds::priority_queue<T, Cmp, __gnu_pbds::rc_binomial_heap_tag>;

template<typename T, typename Cmp = less<T> >
using PairingHeap = __gnu_pbds::priority_queue<T, Cmp, __gnu_pbds::pairing_heap_tag>;

template<typename T, typename Cmp = less<T> >
using ThinHeap = __gnu_pbds::priority_queue<T, Cmp, __gnu_pbds::thin_heap_tag>;

MOLSQ_NAMESPACE_HPP_END

#endif // PQUEUE_HPP