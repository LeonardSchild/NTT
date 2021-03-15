//
// Created by leonard on 13.03.21.
//

#ifndef ENTITY_INTERFACES_H
#define ENTITY_INTERFACES_H

#include <climits>
#include <concepts>
#include <vector>
#include <cstdint>

/** Abstract template class to model a finite field
 * @tparam T variable type for computations
 */
template<typename T>
class FiniteField {

public:
    /**
     *  Constructor
     * @param modulus modulus of finite field
     */
    explicit FiniteField(T modulus) : modulus(modulus), bitwidth(sizeof(T) * CHAR_BIT) {};

    /**
     * Destructor
     */
    virtual ~FiniteField() = default;

    /**
     * Maps elements of the finite field to their Montgomery for modulus 2^\p bitwidth
     * @param input Element to map
     * @return \p input in Montgomery form
     */
    virtual T SwitchTo2N(T input) = 0;

    /**
     * Maps element from their Montgomery form to their representation modulo \p modulus
     * @param input Element in Montgomery form
     * @return \p input modulo \p modulus
     */
    virtual T SwitchFrom2N(T input) = 0;

    /**
     * Multiplies two elements of the finite field given in Montgomery form.
     * @param a An element from the finite field in Montgomery form
     * @param b An element from the finite field in Montgomery form
     * @return The product of \p a and \p b in Montgomery form
     */
    virtual T ModMul(T a, T b) = 0;

    /**
     * Adds two elements of the finite field given in Montgomery form
     * @param a An element from the finite field in Montgomery form
     * @param b An element from the finite field in Montgomery form
     * @return The sum of \p a and \p b
     */
    virtual T ModAdd(T a, T b) = 0;

    /**
     * Subtracts two elements of the finite field given in Montgomery form
     * @param a An element from the finite field in Montgomery form
     * @param b An element from the finite field in Montgomery form
     * @return The difference of \p a and \p b
     */
    virtual T ModSub(T a, T b) = 0;

    /**
     * Raises a element of the finite field given in Montgomery Form to a power.
     * @param a An element from the finite field in Montgomery form
     * @param e The power to which \p a should be raised
     * @return \p a raised to the \p e-th power
     */
    virtual T ModExp(T a, T e) = 0;

    /**
     * Computes a primitive root of unity given the order
     * @param order Order of the primitive root. We assume it divides \p modulus - 1
     * @return A primitive root of unity of given order
     */
    virtual T ComputePrimitiveRootOfUnity(T order) = 0;

    T GetModulus() const { return modulus; }

protected:

    /**
     * Precomputes parameters based on the modulus
     */
    virtual void PreCompute() = 0;

    //! Modulus of finite field
    T modulus;
    //! Negated inverse of \modulus modulo 2^\bitwidth
    T modulus_inverse;
    //! Number of bits used by the type T
    T bitwidth;

};

/**
 * Worker engine to perform Number Theoretic Transforms
 * @tparam T Type to use for computations
 * @tparam FF_class Class modeling the underlying finite field we use for computations
 */
template<typename T, class FF_class> requires std::derived_from<FF_class, FiniteField<T>>
class NttEngine {
public:

    /**
     * Constructor of the engine
     * @param modulus Modulus for the finite field
     * @param N Length of arrays / degree of polynomials we will use
     */
    NttEngine(T modulus, T N) : FF(modulus), N(N) {};

    /**
     * Destructor
     */
    virtual ~NttEngine() = default;

    /**
     * Multiplies each entry by the appropriate power of \p phi and performs a forward NTT
     * @param in_out Input AND output vector of the data to transform. It is assumed that entries are already in montgomery format
     */
    virtual void Forward(std::vector<T>& in_out) = 0;

    /**
     * Multiplies each entry by the appropriate power of the \p phi_inv and performs a forward NTT
     * @param in_out Input AND output vector of the data to transform. It is assumed that entries are already in montgomery format
     */
    virtual void Backward(std::vector<T>& in_out) = 0;

    /**
     * Multiplies two polynomials to which the Forward routine has already been applied to. I.e. we basically compute
     * the Hadamard product
     * @param out Output variable which will contain \p lhs * \p rhs
     * @param lhs Left input of the product
     * @param rhs Right input of the product
     */
    virtual void Multiply(std::vector<T>& out, std::vector<T>& lhs, std::vector<T>& rhs) = 0;

    /**
     * Multiplies two polynomials in "plain" format, i.e. they are not in NTT format yet. We still require that
     * the Montgomery transform has been appliedi to the coefficients.
     * @param out Output variable which will contain \p lhs * \p rhs
     * @param lhs Left input of the product
     * @param rhs Right input of the product
     */
    virtual void TransformAndMultiply(std::vector<T>& out, std::vector<T>& lhs, std::vector<T>& rhs) = 0;

protected:

    /**
     * Enum to specify transform direction
     */
    enum Direction {
        FORWARD = 0,
        BACKWARD = 1
    };

    /**
     * Precomputes coefficients based on the parameters \p N and \p modulus of the constructor.
     */
    virtual void PreCompute() = 0;

    /**
     * Both the forward and backward transform can be defined using a single core routine (plus some additional operations)
     * which is implemented here.
     * @param in_out Input and output of the data to transform.
     * @param d Kind of transform. Forward or Backward.
     */
    virtual void ProtoTransform(std::vector<T>& in_out, Direction d) = 0;

    //! Underlying Finite Field
    FF_class FF;

    //! Transform length and its logarithm base 2
    uint32_t N, logN;
    //! Inverse of transform length w.r.t. the modulus
    T invN;
    //! Pair which contains the 2*N-th root of unity and its inverse w.r.t. the modulus
    T phi[2];
    //! Pair which contains the N-th root of unity and its inverse w.r.t. the modulus
    T omega[2];
};

#endif //ENTITY_INTERFACES_H
