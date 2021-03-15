//
// Created by leonard on 22.02.21.
//

#ifndef ENTITY_UTILS_H
#define ENTITY_UTILS_H

#include <cstdint>
#include <climits>
#include <type_traits>

#include <iostream>
/**
 * Builds a __uint128_t from two uint64_t operands
 * @param a High part of the output
 * @param b Low part of the output
 * @return a * 2^64 + b
 */
inline __uint128_t make_2u64u128(uint64_t a, uint64_t b) {
    return __uint128_t(a) << 64 | b;
}

/**
 * Builds a __uint128_t from four uint32_t operands
 * @return a * 2^96 + b * 2^64 + c * 2^32 + d
 */
inline __uint128_t make_4u32u128(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    return ((__uint128_t(a) << 32) | b) << 64 | ((__uint128_t(c) << 32) | d);
}

/**
 * Prints a __uint128_t
 * @param a The value to be printed
 */
inline void print_u128(__uint128_t &a) {
    uint64_t hi = a >> 64;
    uint64_t lo = (a << 64) >> 64;
    std::cout << "HI: " << hi << " LO: " << lo << std::endl;
}
/**
 * The Clang compiler provides routines two reverse the bits of an operand
 * This templated function is used to pick the appropriate one by type
 * @tparam T Type of the operand
 * @param v The value whose bits shall be reversed
 * @return rev(v)
 */
template<typename T>
T clang_reverse(T v) {
    static_assert(std::is_integral<T>::value);
    if constexpr (sizeof(T) * CHAR_BIT == 8) {
        return __builtin_bitreverse8(v);
    }
    if constexpr (sizeof(T) * CHAR_BIT == 16) {
        return __builtin_bitreverse16(v);
    }
    if constexpr (sizeof(T) * CHAR_BIT == 32) {
        return __builtin_bitreverse32(v);
    } else {
        return __builtin_bitreverse64(v);
    }
}
/**
 * Reverses a variable amount of low order bits from a 32-bit operand and discards the remaining ones.
 * @param v The value to be reversed
 * @param bits The number of low order bits
 * @return rev(v) >> (32 - bits)
 */
inline uint32_t bit_reverse(uint32_t v, uint32_t bits) {
#if defined(__clang__)
    return clang_reverse<uint32_t>(v) >> (sizeof(uint32_t) * CHAR_BIT - bits);
#else
    v = ((v >> 1) & 0x55555555u) | ((v & 0x55555555u) << 1);
    v = ((v >> 2) & 0x33333333u) | ((v & 0x33333333u) << 2);
    v = ((v >> 4) & 0x0f0f0f0fu) | ((v & 0x0f0f0f0fu) << 4);
    v = ((v >> 8) & 0x00ff00ffu) | ((v & 0x00ff00ffu) << 8);
    v = ((v >> 16) & 0xffffu) | ((v & 0xffffu) << 16);
    return v >> (sizeof(uint32_t) * CHAR_BIT - bits);
#endif
}

/**
 * Reverses a variable amount of low order bits from a 64-bit operand and discards the remaining ones.
 * @param v The value to be reversed
 * @param bits The number of low order bits
 * @return rev(v) >> (64 - bits)
 */
inline uint64_t bit_reverse(uint64_t v, uint32_t bits) {
#if defined(__clang__)
    return clang_reverse<uint64_t>(v) >> (sizeof(uint64_t) * CHAR_BIT - bits);
#else
    uint32_t lo = v & 0xffffffff;
    uint32_t hi = v >> 32;
    return ((uint64_t(bit_reverse(lo, sizeof(uint32_t) * CHAR_BIT)) << 32)
    | bit_reverse(hi, sizeof(uint32_t) * CHAR_BIT )) >> (sizeof(uint64_t) * CHAR_BIT - bits);
#endif
}

/**
 * Reverses a variable amount of low order bits from a 128-bit operand and discards the remaining ones.
 * @param v The value to be reversed
 * @param bits The number of low order bits
 * @return rev(v) >> (128 - bits)
 */
inline __uint128_t bit_reverse(__uint128_t v, uint32_t bits) {
    uint64_t lo = v & (0xffffffffffffffff);
    uint64_t hi = v >> 64;
    return ((__uint128_t(bit_reverse(lo, sizeof(uint64_t) * CHAR_BIT)))
            | bit_reverse(hi, sizeof(uint64_t) * CHAR_BIT)) >> (sizeof(__uint128_t) * CHAR_BIT - bits);
}

/**
 * Computes the (ceil of the) base 2 logarithm of a value, either by counting the amount of leading zeroes or
 * by actually computing the logarithm
 * @param v Input value
 * @return ceil(log_2(v))
 */
inline uint32_t mylog2(uint32_t v) {
#if defined(__clang__)
    return sizeof(uint32_t) * CHAR_BIT - __builtin_clz(v);
#elif defined(__GNUG__)
    return sizeof(uint32_t) * CHAR_BIT - __builtin_clz(v);
#else
    return std::log(v) / std::log(2);
#endif
}

#endif //ENTITY_UTILS_H
