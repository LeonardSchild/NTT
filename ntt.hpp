#include "common.cuh"

template<typename T, T modulus, T poly_dim>
class NTT_worker {

	/* See for details
	 * https://en.wikipedia.org/wiki/Root_of_unity_modulo_n#Finding_a_primitive_k-th_root_of_unity_modulo_n
	 */
	public:
		NTT_worker(): modulus(modulus), log_poly_dim(LOG_2<T, poly_dim>::n) {
			// is poly_dim a power of 2 ?
			static_assert((poly_dim & (poly_dim - 1)) == 0);
			// are there square roots for omega ?
			static_assert((modulus % (2 * poly_dim)) == 1);

			// prime modulus carmichael = phi function
			// x^(carmichael(p)/k) = 1 ^ (-k) mod p
			T cm = (modulus - 1) / poly_dim;
			T cm2 = (modulus - 1) / (2 * poly_dim);
			omega = mod_exp<T, modulus>(2, cm);
			phi = mod_exp<T, modulus>(2, cm2);

			// use small fermat for inverse
			// a^(p-1) = 1 mod p => a^(p-2) = a^-1 mod p
			omega_inverse = mod_exp<T, modulus>(omega, modulus - 2);
			phi_inverse = mod_exp<T, modulus>(phi, modulus - 2);
			poly_dim_inverse = mod_exp<T, modulus>(poly_dim, modulus - 2);

			T o = 1, io = 1, p = 1, ip = 1;
			for(int i = 0; i < poly_dim; i++) {
				powers_of_omega[i] = o;
				powers_of_omega_inverse[i] = io;
				powers_of_phi[i] = p;
				powers_of_phi_inverse[i] = ip;

				o = mod_mul<T, modulus>(o, omega);
				io = mod_mul<T, modulus>(io, omega_inverse);
				p = mod_mul<T, modulus>(p, phi);
				ip = mod_mul<T, modulus>(ip, phi_inverse);
			}
		}

		void transform(RingPolynomial<T, poly_dim>& result, RingPolynomial<T, poly_dim>& input) {
			assert(input.state == RAW);
			proto_transform(result, input, powers_of_omega);
			result.state = NTT;
		}

		void inverse_transform(RingPolynomial<T, poly_dim>& result, RingPolynomial<T, poly_dim>& input) {
			assert(input.state == NTT);
			proto_transform(result, input, powers_of_omega_inverse);
			result.state = RAW;

			for(int i = 0; i < poly_dim; i++) {
				result->coefficients[i] = mod_mul<T, modulus>(result->coefficients[i], poly_dim_inverse);
			}
		}

		void multiply(RingPolynomial<T, poly_dim>& result, RingPolynomial<T, poly_dim>& lhs, RingPolynomial<T, poly_dim>& rhs) {

			product_buffer_1 = lhs->coefficients;
			product_buffer_2 = rhs->coefficients;

			for(int i = 0; i < poly_dim; i++) {
				product_buffer_1[i] = mod_mul<T, modulus>(product_buffer_1[i], powers_of_phi[i]);
				product_buffer_2[i] = mod_mul<T, modulus>(product_buffer_2[i], powers_of_phi[i]);
			}

			transform(product_buffer_1, product_buffer_1);
			transform(product_buffer_2, product_buffer_2);

			for(int i = 0; i < poly_dim; i++) {
				result->coefficients[i] = mod_mul<T, modulus>(product_buffer_1[i], product_buffer_2[i]);
			}

			inverse_transform(result->coefficients, result->coefficients);

			for(int i = 0; i < poly_dim; i++) {
				result->coefficients[i] = mod_mul<T, modulus>(result->coefficients[i], powers_of_phi_inverse[i]);
			}
		}

	private:
		void proto_transform(RingPolynomial<T, poly_dim>& result, RingPolynomial<T, poly_dim>& input, std::array<T, poly_dim>& twiddle_factors) {

			transform_buffer.fill(T(0));
			for(uint32_t i = 0; i < poly_dim; i++) {
				transform_buffer[i] = input->coefficients[reverse_bits<T>(i) >> (sizeof(T) * 8 - log_poly_dim)];
			}

			for(uint32_t i = 0; i < log_poly_dim; i++) {
				uint32_t shift = log_poly_dim - i - 1;
				for(uint32_t j = 0; j < poly_dim/2; j++){

					uint32_t p = (j >> shift) << shift;

					T lhs = transform_buffer[2 * j];
					T rhs = mod_mul<T, mod>(transform_buffer[2 * j + 1], twiddle_factors[p]);

					out->coefficients[j] = (lhs + rhs) % mod;
					out->coefficients[j + (poly_dim / 2)] = mod_sub<T, mod>(lhs, rhs);
				}
				transform_buffer = out->coefficients;
			}

		}

		T omega;
		T omega_inverse;
		T phi;
		T phi_inverse;
		T poly_dim_inverse;
		T log_poly_dim;

		std::array<T, poly_dim> transform_buffer;
		std::array<T, poly_dim> product_buffer_1;
		std::array<T, poly_dim> product_buffer_2;

		std::array<T, poly_dim> powers_of_omega;
		std::array<T, poly_dim> powers_of_omega_inverse;
		std::array<T, poly_dim> powers_of_phi;
		std::array<T, poly_dim> powers_of_phi_inverse;
};

