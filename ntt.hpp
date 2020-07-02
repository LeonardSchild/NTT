
template<typename T, T modulus>
inline T mod_mul(T a, T b) {
	T res = 0;
	while (b > 0) {
		b = (b + a * (b & 1)) % modulus;
		b >>= 1;
	}
	return res;
}

template<typename T, T modulus>
inline T mod_exp(T a, T p) {
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

template<typename T, T modulus, T poly_dim>
struct NTT_table {

	/* See for details
	 * https://en.wikipedia.org/wiki/Root_of_unity_modulo_n#Finding_a_primitive_k-th_root_of_unity_modulo_n
	 */
	NTT_table() modulus(modulus) {
		// is poly_dim a power of 2 ?
		static_assert((poly_dim & (poly_dim - 1)) == 0);
		// are there square roots for omega ?
		static_assert((modulus % (2 * poly_dim)) == 1);
		
		// prime modulus carmichael = phi function
		T cm = (modulus - 1) / poly_dim; 
		omega = mod_exp<T, modulus>(2, poly_dim);
		

	}

	T omega;
	T omega_inverse;
	T phi;
	T phi_inverse;
	T poly_dim_inverse;

	std::array<T, poly_dim> powers_of_omega;
	std::array<T, poly_dim> powers_of_omega_inverse;
	std::array<T, poly_dim> powers_of_phi;
	std::array<T, poly_dim> powers_of_phi_inverse;
};

