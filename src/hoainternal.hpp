/**
 * Copyright Tom van Dijk
 */

/*
 * Internals of the eHOA parser (simplehoa.cpp); not part of the interface in
 * simplehoa.h.
 */

#ifndef KNOR_HOAINTERNAL_HPP
#define KNOR_HOAINTERNAL_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "simplehoa.h"

namespace hoa {

inline void* checkedMalloc(size_t bytes)
{
    void* p = malloc(bytes == 0 ? 1 : bytes);
    if (p == nullptr) throw std::bad_alloc();
    return p;
}

inline void* checkedCalloc(size_t count, size_t size)
{
    void* p = calloc(count == 0 ? 1 : count, size);
    if (p == nullptr) throw std::bad_alloc();
    return p;
}

inline char* mallocString(const std::string& s)
{
    char* p = static_cast<char*>(checkedMalloc(s.size() + 1));
    memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

template <typename T>
inline T* mallocArray(const std::vector<T>& v)
{
    if (v.empty()) return nullptr;
    T* p = static_cast<T*>(checkedMalloc(v.size() * sizeof(T)));
    memcpy(p, v.data(), v.size() * sizeof(T));
    return p;
}

// Bump allocator for the parts of an automaton that are freed all at once.
class Arena final {
public:
    Arena() = default;
    Arena(const Arena&) = delete;
    Arena& operator=(const Arena&) = delete;

    template <typename T>
    T* allocate(size_t count)
    {
        static_assert(alignof(T) <= alignof(std::max_align_t), "unsupported alignment");
        const size_t bytes = count * sizeof(T);
        uintptr_t p = alignUp(reinterpret_cast<uintptr_t>(cur), alignof(T));
        const uintptr_t limit = reinterpret_cast<uintptr_t>(end);
        if (cur == nullptr || p > limit || bytes > limit - p) {
            grow(bytes);
            p = reinterpret_cast<uintptr_t>(cur);
        }
        cur = reinterpret_cast<char*>(p + bytes);
        return reinterpret_cast<T*>(p);
    }

    template <typename T>
    T* copy(const T* data, size_t count)
    {
        if (count == 0) return nullptr;
        T* p = allocate<T>(count);
        memcpy(p, data, count * sizeof(T));
        return p;
    }

    char* copyString(const std::string& s)
    {
        char* p = allocate<char>(s.size() + 1);
        memcpy(p, s.c_str(), s.size() + 1);
        return p;
    }

private:
    static constexpr size_t FIRST_BLOCK = size_t(64) << 10;
    static constexpr size_t LARGEST_BLOCK = size_t(16) << 20;

    static uintptr_t alignUp(uintptr_t p, size_t alignment)
    {
        return (p + alignment - 1) & ~uintptr_t(alignment - 1);
    }

    void grow(size_t bytes)
    {
        const size_t size = std::max(bytes, next_block);
        if (next_block < LARGEST_BLOCK) next_block *= 2;
        std::unique_ptr<char[]> block(new char[size]);
        cur = block.get();
        blocks.push_back(std::move(block));
        end = cur + size;
    }

    std::vector<std::unique_ptr<char[]>> blocks;
    char* cur = nullptr;
    char* end = nullptr;
    size_t next_block = FIRST_BLOCK;
};

// Recursively free a malloc'd tree (the representation used outside the arena).
inline void freeTree(BTree* t)
{
    while (t != nullptr) {
        freeTree(t->right);
        free(t->alias);
        BTree* left = t->left;
        free(t);
        t = left;
    }
}

struct TreeDeleter {
    void operator()(BTree* t) const { freeTree(t); }
};
using TreePtr = std::unique_ptr<BTree, TreeDeleter>;

inline TreePtr mallocNode(NodeType type, TreePtr left, TreePtr right, int id)
{
    BTree* n = static_cast<BTree*>(checkedMalloc(sizeof(BTree)));
    n->left = left.release();
    n->right = right.release();
    n->alias = nullptr;
    n->id = id;
    n->type = type;
    return TreePtr(n);
}

inline bool hasString(char* const* list, int count, const char* s)
{
    for (int i = 0; i < count; i++) {
        if (strcmp(list[i], s) == 0) return true;
    }
    return false;
}

// Is t exactly Inf(p) (if good) or Fin(p) (otherwise)?
inline bool isParityTerm(const BTree* t, bool good, int p)
{
    if (t == nullptr || t->type != (good ? NT_INF : NT_FIN)) return false;
    const BTree* set = t->left;
    return set != nullptr && set->type == NT_SET && set->id == p;
}

// Match the canonical parity condition from priority p onwards, e.g.
// "parity max even 4": Fin(3) & (Inf(2) | (Fin(1) & Inf(0))). Unless strict,
// this is bug-compatible with the original hoa-tools check of isParityGFG:
// in a disjunction, one matching side suffices.
inline bool isParityCondition(const BTree* t, int noSets, bool max, short good, int p, bool strict)
{
    if (t == nullptr) return false;
    const bool isGood = p % 2 == good;
    const bool last = max ? p == 0 : p == noSets - 1;
    if (last) return isParityTerm(t, isGood, p);
    if (t->type != (isGood ? NT_OR : NT_AND)) return false;
    const bool head = isParityTerm(t->left, isGood, p);
    const bool tail = isParityCondition(t->right, noSets, max, good, max ? p - 1 : p + 1, strict);
    return isGood && !strict ? (head || tail) : (head && tail);
}

// Read "parity (max|min) (even|odd)" from acc-name: 0 if found, otherwise
// the error code of isParityGFG.
inline int parityName(const HoaData* data, bool* isMaxParity, short* resGoodPriority)
{
    if (data->accNameID == nullptr || strcmp(data->accNameID, "parity") != 0) return 100;
    bool foundOrder = false;
    bool foundGood = false;
    for (int i = 0; i < data->noANPs; i++) {
        const char* p = data->accNameParameters[i];
        if (strcmp(p, "max") == 0 || strcmp(p, "min") == 0) {
            *isMaxParity = p[1] == 'a';
            foundOrder = true;
        } else if (strcmp(p, "even") == 0 || strcmp(p, "odd") == 0) {
            *resGoodPriority = p[0] == 'e' ? 0 : 1;
            foundGood = true;
        }
    }
    if (!foundOrder) return 101;
    if (!foundGood) return 102;
    return 0;
}

}  // namespace hoa

#endif
