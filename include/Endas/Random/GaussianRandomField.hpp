/**
 * @file GaussianRandomField.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_RANDOM_GAUSSIAN_RANDOM_FIELD_HPP__
#define __ENDAS_RANDOM_GAUSSIAN_RANDOM_FIELD_HPP__


#include <Endas/Spatial/Variogram.hpp>

#include <memory>


namespace endas
{

/** 
 * @addtogroup random
 * @{ 
 */


/**
 * Stationary Gaussian Random Field (GRF) generator in two dimensions.
 * 
 * The generator uses the circulant embedding method so that the usual sampling by Cholesky
 * decomposition of full-rank covariance matrices can be replaced by a Fourier transform.
 * This allows much larger random fields to be generated than what is possible with the Cholesky
 * approach (see the MultivariateRandomNormal class).
 * 
 * However, the embedding may fail for some covariance functions, leading to negative eigenvalues.
 * One solution is to increase the size of the matrix in which the covariance is embedded until
 * eigenvalues are positive or zero. In this implementation, the enlargement is done if any
 * eigenvalue is smaller than ``eigMin`` and the matrix is of reasonable size. Since any (slightly)
 * negative eigenvalues will at the end be pruned, there may be artifacts in the output if large
 * (in the negative sense) value of ``eigMin`` is allowed. Setting ``eigMin`` to zero, on the other
 * hand, is seldom needed and will likely lead to significant decrease in performance.
 * 
 * Because the implementation relies on pre-compued fft, the size of the random field to generate 
 * must be given to the constructor 
 */
class ENDAS_DLL GaussianRandomField
{
public:

    /** 
     * GaussianRandomField constructor.
     * 
     * @param ny       Size of the generated random field in y dimension
     * @param nx       Size of the generated random field in x dimension
     * @param covFn    Covariance function
     * @param eigMin   Minimum allowed eigenvalues of the circulant embedding of the covariance.
     */
    GaussianRandomField(int nx, int ny, const IsotropicCovarianceFn& covFn, double eigMin = -1e-3);

    ~GaussianRandomField();

    /** 
     * Fills the given array with samples from the distribution.
     */
    void operator()(Ref<Array2d> out) const;


    /** 
     * Fills the given array with samples from the distribution.
     * 
     * This overload generates two realizations of the random field at once. Because two fields 
     * are always generated simultaneuosly, this is more efficient when more than one realization
     * is needed.
     */
    void operator()(Ref<Array2d> out, Ref<Array2d> out2) const;



private:
    struct Data;
    std::unique_ptr<Data> mData;
};


/** @} */

}

#endif