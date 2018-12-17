/*
 * Copyright (C) 2015-2018, Nils Moehrle
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ACC_MATH_HEADER
#define ACC_MATH_HEADER

#include <cstdint>

#include "defines.h"

ACC_NAMESPACE_BEGIN

template <template<typename, int> class VecType>
std::uint64_t z_order_index(VecType<std::uint32_t, 3> position) {
    /* Derived from an implementation by Jens Schneider */
    std::uint64_t ret = 0;
    std::uint64_t x = position[0] & 0x3FFFFFull; // 22 bit range
    x = (x * 0x100000001ull) & 0x003F00000000FFFFull;
    x = (x * 0x10001ull) & 0x003F0000FF0000FFull;
    x = (x * 0x101ull) & 0x300F00F00F00F00Full;
    x = (x * 0x11ull) & 0x30C30C30C30C30C3ull;
    ret = (x * 0x5ull) & 0x9249249249249249ull;

    std::uint64_t y = position[1] & 0x1FFFFFull; // 21 bit range
    y = (y * 0x100000001ull) & 0x003F00000000FFFFull;
    y = (y * 0x10001ull) & 0x003F0000FF0000FFull;
    y = (y * 0x101ull) & 0x300F00F00F00F00Full;
    y = (y * 0x11ull) & 0x30C30C30C30C30C3ull;
    y = (y * 0x5ull) & 0x9249249249249249ull;
    ret |= (y << 1ull);

    std::uint64_t z = position[2] & 0x1FFFFFull; // 21 bit range
    z = (z * 0x100000001ull) & 0x003F00000000FFFFull;
    z = (z * 0x10001ull) & 0x003F0000FF0000FFull;
    z = (z * 0x101ull) & 0x300F00F00F00F00Full;
    z = (z * 0x11ull) & 0x30C30C30C30C30C3ull;
    z = (z * 0x5ull) & 0x9249249249249249ull;
    ret |= (z << 2ull);

    return ret;
}

template <template<typename, int> class VecType>
VecType<std::uint32_t, 3> z_order_position(std::uint64_t index) {
    /* Derived from an implementation by Jens Schneider */
    VecType<std::uint64_t, 3> tmp;
    VecType<std::uint32_t, 3> ret;
    tmp[0] = index & 0x9249249249249249ull;
    tmp[0] = (tmp[0] | (tmp[0] >> 2ull)) & 0x30C30C30C30C30C3ull;
    tmp[0] = (tmp[0] | (tmp[0] >> 4ull)) & 0x300F00F00F00F00Full;
    tmp[0] = (tmp[0] | (tmp[0] >> 8ull)) & 0x003F0000FF0000FFull;
    tmp[0] = (tmp[0] | (tmp[0] >> 16ull)) & 0x003F00000000FFFFull;
    tmp[0] = (tmp[0] | (tmp[0] >> 32ull)) & 0x3FFFFFull;
    ret[0] = std::uint32_t(tmp[0]);

    tmp[1] = (index >> 1ull) & 0x9249249249249249ull;
    tmp[1] = (tmp[1] | (tmp[1] >> 2ull)) & 0x30C30C30C30C30C3ull;
    tmp[1] = (tmp[1] | (tmp[1] >> 4ull)) & 0x300F00F00F00F00Full;
    tmp[1] = (tmp[1] | (tmp[1] >> 8ull)) & 0x003F0000FF0000FFull;
    tmp[1] = (tmp[1] | (tmp[1] >> 16ull)) & 0x003F00000000FFFFull;
    tmp[1] = (tmp[1] | (tmp[1] >> 32ull)) & 0x3FFFFFull;
    ret[1] = std::uint32_t(tmp[1]);

    tmp[2] = (index >> 2ull) & 0x9249249249249249ull;
    tmp[2] = (tmp[2] | (tmp[2] >> 2ull)) & 0x30C30C30C30C30C3ull;
    tmp[2] = (tmp[2] | (tmp[2] >> 4ull)) & 0x300F00F00F00F00Full;
    tmp[2] = (tmp[2] | (tmp[2] >> 8ull)) & 0x003F0000FF0000FFull;
    tmp[2] = (tmp[2] | (tmp[2] >> 16ull)) & 0x003F00000000FFFFull;
    tmp[2] = (tmp[2] | (tmp[2] >> 32ull)) & 0x3FFFFFull;
    ret[2] = std::uint32_t(tmp[2]);

    return ret;
}

ACC_NAMESPACE_END

#endif /* ACC_MATH_HEADER */
