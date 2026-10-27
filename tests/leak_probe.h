#ifndef LEAK_PROBE_H
#define LEAK_PROBE_H

/*
 * leak_probe.h -- portable leak detector for the LS-align test harness.
 *
 * ASan/UBSan are NOT available on the MSYS2/MinGW toolchain used here
 * (the libasan runtime is missing; linking with -fsanitize=address fails).
 * Instead we replace the global operator new/delete inside the *test*
 * translation unit only -- the shipping source is never modified -- and count
 * live allocation events.  A scope that leaks shows a positive delta of
 * liveAllocs().
 *
 * This is a coarse instrument: it counts allocation events, not bytes, and it
 * cannot see C-level malloc()/free() inside the library.  It is precise enough
 * for the allocation leaks in LS-align (AlignInfo, Dijkstra's flag array,
 * std::vector/std::string held by leaked objects).
 *
 * NOTE: the toolchain's static libstdc++ is built for C++17 and may route some
 * allocations through the aligned operator new/delete overloads, so those are
 * provided as well.  Compile this TU as C++17 to match.
 */

#include <cstdlib>
#include <cstddef>
#include <new>

namespace leakprobe {
inline long& live() { static long n = 0; return n; }
inline long& mismatch() { static long n = 0; return n; }

const unsigned SCALAR_MAGIC = 0x5343414cu;
const unsigned ARRAY_MAGIC  = 0x41525241u;

struct Header { unsigned magic; unsigned pad[3]; };
}

inline long liveAllocs() { return leakprobe::live(); }
inline long mismatchCount() { return leakprobe::mismatch(); }

/* Replacement allocator: route everything through malloc/free so the count is
 * symmetric across all overloads used by this single-TU program. */
inline void* lp_alloc(std::size_t n, unsigned tag) {
    leakprobe::Header* h = static_cast<leakprobe::Header*>(
        std::malloc(sizeof(leakprobe::Header) + (n ? n : 1)));
    if (!h) return 0;
    h->magic = tag;
    ++leakprobe::live();
    return h + 1;
}
inline void lp_free(void* p, unsigned expected) {
    if (!p) return;
    leakprobe::Header* h = static_cast<leakprobe::Header*>(p) - 1;
    if (h->magic != expected) ++leakprobe::mismatch();
    --leakprobe::live();
    std::free(h);
}

void* operator new(std::size_t n) { void* p = lp_alloc(n, leakprobe::SCALAR_MAGIC); if (!p) throw std::bad_alloc(); return p; }
void* operator new[](std::size_t n) { void* p = lp_alloc(n, leakprobe::ARRAY_MAGIC); if (!p) throw std::bad_alloc(); return p; }
void* operator new(std::size_t n, const std::nothrow_t&) noexcept { return lp_alloc(n, leakprobe::SCALAR_MAGIC); }
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept { return lp_alloc(n, leakprobe::ARRAY_MAGIC); }

void operator delete(void* p) noexcept { lp_free(p, leakprobe::SCALAR_MAGIC); }
void operator delete[](void* p) noexcept { lp_free(p, leakprobe::ARRAY_MAGIC); }
void operator delete(void* p, const std::nothrow_t&) noexcept { lp_free(p, leakprobe::SCALAR_MAGIC); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { lp_free(p, leakprobe::ARRAY_MAGIC); }
void operator delete(void* p, std::size_t) noexcept { lp_free(p, leakprobe::SCALAR_MAGIC); }
void operator delete[](void* p, std::size_t) noexcept { lp_free(p, leakprobe::ARRAY_MAGIC); }

#if defined(__cpp_aligned_new)
void* operator new(std::size_t n, std::align_val_t) { void* p = lp_alloc(n, leakprobe::SCALAR_MAGIC); if (!p) throw std::bad_alloc(); return p; }
void* operator new[](std::size_t n, std::align_val_t) { void* p = lp_alloc(n, leakprobe::ARRAY_MAGIC); if (!p) throw std::bad_alloc(); return p; }
void operator delete(void* p, std::align_val_t) noexcept { lp_free(p, leakprobe::SCALAR_MAGIC); }
void operator delete[](void* p, std::align_val_t) noexcept { lp_free(p, leakprobe::ARRAY_MAGIC); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { lp_free(p, leakprobe::SCALAR_MAGIC); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { lp_free(p, leakprobe::ARRAY_MAGIC); }
#endif

/* RAII helper: measure the change in live allocations over a scope. */
struct LeakDelta {
    long base;
    LeakDelta() : base(leakprobe::live()) {}
    long delta() const { return leakprobe::live() - base; }
};

struct MismatchDelta {
    long base;
    MismatchDelta() : base(leakprobe::mismatch()) {}
    long delta() const { return leakprobe::mismatch() - base; }
};

#endif /* LEAK_PROBE_H */
