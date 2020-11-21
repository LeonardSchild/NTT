#ifndef NTT_CUH
#define NTT_CUH

#include "common.cuh"
#include "hontt.cuh"
#include "parameters.cuh"
#include "utils.cuh"


template <typename T, T mod>
__global__ void mult_vector(T* out, T *lhs, T *rhs, size_t offset = 0) {

	uint32_t in_idx = threadIdx.x;
	uint32_t el_idx = blockIdx.x * blockDim.x;

	(out + el_idx)[in_idx] = mod_mul<T, mod>((lhs + el_idx)[in_idx], (rhs + el_idx)[in_idx]);

}

template <typename T, T mod, T dim>
__global__ void mult_phi(NTT_table<T,dim>* table, T* out, T* in) {

	uint32_t in_idx = threadIdx.x;
	uint32_t el_idx = blockIdx.x * blockDim.x;

	(out + el_idx)[in_idx] = mod_mul<T, mod>((in + el_idx)[in_idx], table->powers_of_phi[in_idx]);
}

template <typename T, T mod, T dim>
__global__ void mult_inv_phi(NTT_table<T,dim>* table, T* out, T* in) {

	uint32_t in_idx = threadIdx.x;
	uint32_t el_idx = blockIdx.x * blockDim.x;

	(out + el_idx)[in_idx] = mod_mul<T, mod>((in + el_idx)[in_idx], table->powers_of_iphi[in_idx]);
}

template <typename T, T mod, T dim, T dim_inverse_mod>
__global__ void proto_transform(T* out_poly, T* in_poly, T* twiddle_factors) {
	__shared__ T buffer[dim];

	constexpr T log_dim = LOG_2<T, dim>::n;

	uint32_t idx = threadIdx.x;
	T* in = total_in + 2 * blockIdx.x * blockDim.x;
	T* out = total_out + 2 * blockIdx.x * blockDim.x;

	buffer[2 * idx] = in_poly[reverse_bits(2 * idx) >> (sizeof(T) * 8 - log_dim)];
	buffer[2 * idx + 1] = in_poly[__brev(2 * idx + 1) >> (sizeof(T) * 8 - log_dim)];

	__syncthreads();

	for(uint32_t i = 0; i < log_dim; i++) {

		uint32_t pij = (idx >> (log_dim - i - 1)) << (log_dim - i - 1);
		T lhs = buffer[2 * idx];
		T rhs = mod_mul<T, mod>(buffer[2 * idx + 1], twiddle_factors[pij]);

		__syncthreads();

		buffer[idx] = (lhs + rhs) % mod; // assume mod < 2^31 or 2^63
		buffer[idx + (dim/2)] = mod_sub<T, mod>(lhs, rhs);

		__syncthreads();
	}

	if constexpr(dim_inverse_mod == 0) {
		out[2*idx] = buffer[2*idx];
		out[2*idx + 1] = buffer[2*idx + 1];
	} else {
		out_poly[2 * idx] = mod_mul<T, mod>(dim_inverse_mod, buffer[2 * idx]);
		out_poly[2 * idx] = mod_mul<T, mod>(dim_inverse_mod, buffer[2 * idx]);
	}
}

template <typename T, T mod, T dim, T log_dim>
__global__ void gpu_NTT(NTT_table<T,dim>* table, T* total_out, T* total_in) {

	__shared__ T buffer[dim];

	uint32_t idx = threadIdx.x;
	T* in = total_in + 2 * blockIdx.x * blockDim.x;
	T* out = total_out + 2 * blockIdx.x * blockDim.x;

	buffer[2 * idx] = in[__brev(2 * idx) >> (32 - log_dim)];
	buffer[2 * idx + 1] = in[__brev(2 * idx + 1) >> (32 - log_dim)];

	__syncthreads();

	for(uint32_t i = 0; i < log_dim; i++) {

		uint32_t pij = (idx >> (log_dim - i - 1)) << (log_dim - i - 1);
		T lhs = buffer[2 * idx];
		T rhs = mod_mul<T, mod>(buffer[2 * idx + 1], table->powers_of_omega[pij]);

		__syncthreads();

		buffer[idx] = (lhs + rhs) % mod; // assume mod < 2^31 or 2^63
		buffer[idx + (dim/2)] = mod_sub<T, mod>(lhs, rhs);

		__syncthreads();
	}

	out[2*idx] = buffer[2*idx];
	out[2*idx + 1] = buffer[2*idx + 1];
}

template <typename T, T mod, T inv_N, uint32_t dim, uint32_t log_dim>
__global__ void gpu_INTT(* table, T* total_out, T* total_in) {

	__shared__ T buffer[dim];

	uint32_t idx = threadIdx.x;
	T* in = total_in + 2 * blockIdx.x * blockDim.x;
	T* out = total_out + 2 * blockIdx.x * blockDim.x;

	buffer[2 * idx] = in[__brev(2 * idx) >> (32 - log_dim)];
	buffer[2 * idx + 1] = in[__brev(2 * idx + 1) >> (32 - log_dim)];

	__syncthreads();

	for(uint32_t i = 0; i < log_dim; i++) {

		uint32_t pij = (idx >> (log_dim - i - 1)) << (log_dim - i - 1);
		T lhs = buffer[2 * idx];
		T rhs = mod_mul<T, mod>(buffer[2 * idx + 1], table->powers_of_iomega[pij]);

		__syncthreads();

		buffer[idx] = (lhs + rhs) % mod;
		buffer[idx + (dim/2)] = mod_sub<T, mod>(lhs, rhs);

		__syncthreads();
	}

	out[2*idx] = mod_mul<T, mod>(buffer[2*idx],inv_N);
	out[2*idx + 1] = mod_mul<T, mod>(buffer[2*idx + 1],inv_N);
}


#endif
