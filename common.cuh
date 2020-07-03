#ifndef __CUDA_ARCH__
#define __device__
#define __host__
#endif

enum FORMAT {
	RAW,
	NTT
};

template<typename T, T poly_dim>
struct RingPolynomial {
	std::array<T, poly_dim> coefficients;
	FORMAT state = NTT;
};

template<typename T, T modulus>
__device__ __host__ inline T mod_sub(T a, T b) {
	return (a >= b) ? a - b : (modulus - b) + a;
}

// maybe inline
template<typename T, T modulus>
__device__ __host__ T mod_mul(T a, T b) {
	T res = 0;
	while (b > 0) {
		b = (b + a * (b & 1)) % modulus;
		b >>= 1;
	}
	return res;
}

//maybe inline
template<typename T, T modulus>
__device__ __host__ T mod_exp(T a, T p) {
	T res = 1
		T accu = a;
	while (p > 0) {
		if (p & 1) {
			res = mod_mul<T, modulus>(res, accu);
		}
		accu = mod_mul<T, modulus>(accu, accu);
		p >>= 1;
	}
}

// https://stackoverflow.com/questions/63776/bit-reversal-of-an-integer-ignoring-integer-size-and-endianness
template<typename T>
__device__ __host__ T reverse_bits_no_specials(T in) {
	T nrev = 0, i, bit1, bit2;

	int count;
	constexpr int type_bits = sizeof(T) * 8;

	for(i = 0; i < type_bits; i += 2)
	{
		/*In each iteration, we  swap one bit on the 'right half'
		  of the number with another on the left half*/

		count =  type_bits - i - 1;  /*this is used to find how many positions
					       to the left (and right) we gotta move
					       the bits in this iteration*/

		bit1 = n & (1<<(i/2)); /*Extract 'right half' bit*/
		bit1 <<= count;         /*Shift it to where it belongs*/

		bit2 = n & 1<<((i/2) + count);  /*Find the 'left half' bit*/
		bit2 >>= count;         /*Place that bit in bit1's original position*/

		nrev |= bit1;   /*Now add the bits to the reversal result*/
		nrev |= bit2;
	}
	return nrev;
}
}

template<typename T>
__device__ __host__ T reverse_bits(T in) {

#if defined(__CUDA_ARCH__)

	if constexpr(sizeof(T) == 4)
		return __brev(in);
	else if constexpr(sizeof(T) == 8)
		return __brevll(in);
	else
		return reverse_bits_no_specials<T>(in);

#elif defined(__CLANG__)
	if constexpr(sizeof(T) == 1)
		return __builtin_bitreverse8(in);
	else if constexpr(sizeof(T) == 2)
		return __builtin_bitreverse16(in);
	else if constexpr(sizeof(T) == 4)
		return __builtin_bitreverse32(in);
	else if constexpr(sizeof(T) == 8)
		return __builtin_bitreverse64(in);
	// maybe for future ints ? __uint128 is a thing after all...
	else
		return reverse_bits_no_specials<T>(in);
#else
	return reverse_bits_no_specials<T>(in);
#endif

}
