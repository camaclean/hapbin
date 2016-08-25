/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014  Colin MacLean <s0838159@sms.ed.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef HAPBIN_HPP
#define HAPBIN_HPP

#include <type_traits>
#include <limits>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include "config.h"

#if defined(__GNUG__)
typedef unsigned long long v8ul  __attribute__((vector_size(64)));
typedef unsigned int       v16ui __attribute__((vector_size(64)));
typedef unsigned long long v4ul  __attribute__((vector_size(32)));
typedef unsigned int       v8ui  __attribute__((vector_size(32)));
typedef unsigned long long v2ul  __attribute__((vector_size(16)));
typedef unsigned int       v4ui  __attribute__((vector_size(16)));
#endif

#ifdef VECLEN
#define VEC VECLEN
#else
#if defined(__AVX512F__) && defined(__GNUG__)
#define VEC 8
#include <zmmintrin.h>
#include <immintrin.h>
#elif defined(__AVX__) && defined(__GNUG__)
#define VEC 4
#elif defined(__SSE2__) && defined(__GNUG__)
#define VEC 2
#else
#define VEC 1
#endif
#endif

#if VEC==8
#define EQUAL(v1, v2) (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3] && v1[4] == v2[4] && v1[5] == v2[5] && v1[6] == v2[6] && v1[7] == v2[7])
#define ZERO ((v8ul){0ULL, 0ULL, 0ULL, 0ULL,0ULL, 0ULL, 0ULL, 0ULL})
#define BITSET_T_MAX ((v8ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__})
#define POPCOUNT popcount8
#elif VEC==4
#define EQUAL(v1, v2) (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3])
#define ZERO ((v4ul){0ULL, 0ULL, 0ULL, 0ULL}) 
#define BITSET_T_MAX ((v4ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__})
#define POPCOUNT popcount4
#elif VEC==2
#define EQUAL(v1, v2) (v1[0] == v2[0] && v1[1] == v2[1])
#define ZERO ((v2ul){0ULL, 0ULL}) 
#define BITSET_T_MAX ((v2ul){__UINT64_MAX__,__UINT64_MAX__})
#define POPCOUNT popcount2
#else
#define EQUAL(v1, v2) (v1 == v2)
#define ZERO 0ULL
#define BITSET_T_MAX __UINT64_MAX__
#define POPCOUNT popcount1
#endif

#if defined(__MINGW32__) && !defined(_ISOC11_SOURCE)
void* aligned_alloc(size_t alignment, size_t size);
void aligned_free(void* ptr);
#elif (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !defined(_ISOC11_SOURCE)
void* aligned_alloc(size_t alignment, size_t size);
#define aligned_free free
#else
#define aligned_free free
#endif

template<typename T, std::size_t N>
constexpr std::size_t ctcBitsetSize()
{
    return (N/(sizeof(T)*8)) + ((N % (sizeof(T)*8) != 0) ? 1 : 0);
}

template<typename T>
std::size_t bitsetSize(std::size_t length)
{
    std::size_t bufferSize = length/(sizeof(T)*8);
    if (length % (sizeof(T)*8) != 0)
        bufferSize += 1;
    return bufferSize;
}

template<typename T, std::size_t N>
constexpr int ctcBitsetMaskShift()
{
    return (N == 0) ? 0 : (sizeof(T)*8 - (N % (sizeof(T)*8)));
}

template<typename T, std::size_t N>
constexpr T ctcBitsetMask()
{
    return (std::numeric_limits<T>::max() >> ctcBitsetMaskShift<T,N>());
}

template<typename T>
inline T bitsetMask(int length)
{
    return (std::numeric_limits<T>::max() >> (sizeof(T)*8-(length % (sizeof(T)*8))));
}

#if VEC==8
inline v8ul bitsetMask8(int length)
{
    int len = length % 512;
    if (len < 64)
        return (v8ul){bitsetMask<unsigned long long>(len),0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL};
    else if (len < 128)
        return (v8ul){__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL,0ULL,0ULL,0ULL,0ULL,0ULL};
    else if (len < 192)
        return (v8ul){__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL,0ULL,0ULL,0ULL,0ULL};
    else if (len < 256)
        return (v8ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL,0ULL,0ULL,0ULL};
    else if (len < 320)
        return (v8ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL,0ULL,0ULL};
    else if (len < 384)
        return (v8ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL,0ULL};
    else if (len < 448)
        return (v8ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL};
    else
        return (v8ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len)};
}
#elif VEC == 4
inline v4ul bitsetMask4(int length)
{
    int len = length % 256;
    if (len < 64)
        return (v4ul){bitsetMask<unsigned long long>(len),0ULL,0ULL,0ULL};
    else if (len < 128)
        return (v4ul){__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL,0ULL};
    else if (len < 192)
        return (v4ul){__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len),0ULL};
    else
        return (v4ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,bitsetMask<unsigned long long>(len)};
}
#elif VEC == 2
inline v2ul bitsetMask2(int length)
{
    int len = length % 128;
    if (len < 64)
        return (v2ul){bitsetMask<unsigned long long>(len),0ULL};
    else
        return (v2ul){__UINT64_MAX__,bitsetMask<unsigned long long>(len)};
}
#endif

template<typename T, typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value>::type* = nullptr>
void convert(const char* line, T* buffer, std::size_t maxLength = 0)
{
    const std::size_t bits = sizeof(T)*8;
    
    std::size_t position = 0;
    for (std::size_t i = 0; line[i] != '\0'; ++i) {
        if (maxLength > 0 && position > maxLength)
            break;
        switch(line[i]) {
            case '0':
                position++;
                break;
            case '1':
                
                buffer[position/bits] |= ((T) (1ULL << (position % bits)));
                position++;
                break;
            case '\n':
            case ' ':
                break;
            default:
                std::cout << "ERROR: Invalid character in ASCII haplotype map: " << line[i] << std::endl;
                throw std::runtime_error("ERROR: Not a valid ASCII haplotype map.");
                break;
        }
    
    }
}

struct Stats 
{
    Stats() : mean(0.0), stddev(0.0) {}
    double mean;
    double stddev;
};

double binom_2(double n);
double nearest(double target, double number);
Stats stats(const std::vector<double>& list);
std::vector<std::string> splitString(const std::string input, char delim);

inline int popcount1(unsigned long long val)
{
#ifdef USE_POPCNT
    return __builtin_popcountll(val);
#else
    val =  val - ((val >> 1) & 0x5555555555555555ULL);
    val = (val & 0x3333333333333333ULL) + ((val >> 2) & 0x3333333333333333ULL);
    val = (val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
    val = (val * 0x0101010101010101ULL) >> 56;
    return val;
#endif
}

inline int popcount8(v8ul val)
{
#ifdef USE_POPCNT
    int ret = __builtin_popcountll(val[0]);
    ret += __builtin_popcountll(val[1]);
    ret += __builtin_popcountll(val[2]);
    ret += __builtin_popcountll(val[3]);
    ret += __builtin_popcountll(val[4]);
    ret += __builtin_popcountll(val[5]);
    ret += __builtin_popcountll(val[6]);
    ret += __builtin_popcountll(val[7]);
    return ret;
#else
    //v8ul  m3    = (v8ul){0x3333333333333333,0x3333333333333333,0x3333333333333333,0x3333333333333333,0x3333333333333333,0x3333333333333333,0x3333333333333333,0x3333333333333333};
    //v8ul  m5    = (v8ul){0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555};
    //v8ul  mF    = (v8ul){0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F};
    //v8ul  mF2   = (v8ul){0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF};
    //v8ul  mF4   = (v8ul){0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF};
    //v8ul  mF8   = (v8ul){0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF};
    //val = (m5  & val) + (m5  & (val >> 1 ));
    //val = (m3  & val) + (m3  & (val >> 2 ));
    //val = (mF  & val) + (mF  & (val >> 4 ));
    //val = (mF2 & val) + (mF2 & (val >> 8 ));
    //val = (mF4 & val) + (mF4 & (val >> 16));
    // //val = (mF8 & val) + (mF8 & (val >> 32));
    //return _mm512_reduce_add_epi32((_m512i) val);
    //v16ui tmp = (v16ui) val;
    //return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7] + tmp[8] + tmp[9] + tmp[10] + tmp[11] + tmp[12] + tmp[13] + tmp[14] + tmp[15];
    val =  val - ((val >> 1) & 0x5555555555555555ULL);
    val = (val & 0x3333333333333333ULL) + ((val >> 2) & 0x3333333333333333ULL);
    val = (val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
    val = (val * 0x0101010101010101ULL) >> 56;
    return val[0] + val[1] + val[2] + val[3] + val[4] + val[5] + val[6] + val[7];
#endif
}

inline int popcount4(v4ul val)
{
#ifdef USE_POPCNT
    int ret = __builtin_popcountll(val[0]);
    ret += __builtin_popcountll(val[1]);
    ret += __builtin_popcountll(val[2]);
    ret += __builtin_popcountll(val[3]);
    return ret;
#else
    //v4ul  m3    = (v4ul){0x3333333333333333,0x3333333333333333,0x3333333333333333,0x3333333333333333};
    //v4ul  m5    = (v4ul){0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555};
    //v4ul  mF    = (v4ul){0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F};
    //v4ul  mF2   = (v4ul){0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF,0x00FF00FF00FF00FF};
    //v4ul  mF4   = (v4ul){0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF,0x0000FFFF0000FFFF};
    //v4ul  mF8   = (v4ul){0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF,0x00000000FFFFFFFF};
    //val = (m5  & val) + (m5  & (val >> 1 ));
    //val = (m3  & val) + (m3  & (val >> 2 ));
    //val = (mF  & val) + (mF  & (val >> 4 ));
    //val = (mF2 & val) + (mF2 & (val >> 8 ));
    //val = (mF4 & val) + (mF4 & (val >> 16));
    //  //val = (mF8 & val) + (mF8 & (val >> 32));
    //v8ui tmp = (v8ui) val;
    //return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
    val =  val - ((val >> 1) & 0x5555555555555555ULL);
    val = (val & 0x3333333333333333ULL) + ((val >> 2) & 0x3333333333333333ULL);
    val = (val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
    val = (val * 0x0101010101010101ULL) >> 56;
    return val[0] + val[1] + val[2] + val[3];
#endif
}

inline int popcount2(v2ul val)
{
#ifdef USE_POPCNT
    int ret = __builtin_popcountll(val[0]);
    ret += __builtin_popcountll(val[1]);
    return ret;
#else
    //v2ul  m3    = (v2ul){0x3333333333333333,0x3333333333333333};
    //v2ul  m5    = (v2ul){0x5555555555555555,0x5555555555555555};
    //v2ul  mF    = (v2ul){0x0F0F0F0F0F0F0F0F,0x0F0F0F0F0F0F0F0F};
    //v2ul  mF2   = (v2ul){0x00FF00FF00FF00FF,0x00FF00FF00FF00FF};
    //v2ul  mF4   = (v2ul){0x0000FFFF0000FFFF,0x0000FFFF0000FFFF};
    //v2ul  mF8   = (v2ul){0x00000000FFFFFFFF,0x00000000FFFFFFFF};
    //val = (m5  & val) + (m5  & (val >> 1 ));
    //val = (m3  & val) + (m3  & (val >> 2 ));
    //val = (mF  & val) + (mF  & (val >> 4 ));
    //val = (mF2 & val) + (mF2 & (val >> 8 ));
    //  //val = (mF4 & val) + (mF4 & (val >> 16));
    val =  val - ((val >> 1) & 0x5555555555555555ULL);
    val = (val & 0x3333333333333333ULL) + ((val >> 2) & 0x3333333333333333ULL);
    val = (val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
    val = (val * 0x0101010101010101ULL) >> 56;
    return val[0] + val[1];
    //val = (mF8 & val) + (mF8 & (val >> 32));
    //v4ui tmp = (v4ui) val;
    //return tmp[0] + tmp[1] + tmp[2] + tmp[3];
#endif
}

#endif // HAPBIN_HPP
