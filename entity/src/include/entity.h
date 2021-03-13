//
// Created by leonard on 13.03.21.
//

#ifndef ENTITY_ENTITY_H
#define ENTITY_ENTITY_H

#include "basic_finite_fields.h"
#include <ntt.h>

using FiniteField32 = BasicFiniteField<uint32_t, uint64_t>;
using NttEngine32 = BasicNttEngine<uint32_t, uint64_t, true>;

#ifdef __SIZEOF_INT128__

using FiniteField64 = BasicFiniteField<uint64_t, __uint128_t>;
using NttEngine64 = BasicNttEngine<uint64_t, __uint128_t, true>;

#endif

#endif //ENTITY_ENTITY_H
