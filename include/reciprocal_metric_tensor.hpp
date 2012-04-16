#ifndef _RECIPROCAL_METRIC_TENSOR_HPP_INCLUDED_DFOSIJEOIUSDFLKJALDFOI4ELKJDLKJFSAINTOSOIJTLKASOFIDJ3LKIJHSFALKJFDL
#define _RECIPROCAL_METRIC_TENSOR_HPP_INCLUDED_DFOSIJEOIUSDFLKJALDFOI4ELKJDLKJFSAINTOSOIJTLKASOFIDJ3LKIJHSFALKJFDL

#include <matrix.hpp>

#include <cmath>

namespace feng
{
namespace reciprocal_metric_tensor_private
{
    template<typename T>
    T 
    f( const T cosalpha, const T cosbeta, const T cosgamma )
    {
        return cosalpha * cosbeta - cosgamma;
    }
}//namespace reciprocal_metric_tensor_private

template<typename T>
const matrix<T>
reciprocal_metric_tensor( const T a, const T b, const T c, const T alpha, const T beta, const T gamma )
{
    using reciprocal_metric_tensor_private::f;

    const T cosalpha = std::cos(alpha);
    const T cosbeta  = std::cos(beta);
    const T cosgamma = std::cos(gamma);
    const T ccc = cosalpha * cosbeta * cosgamma;
    const T abc = a*b*c;
    const T aa = a*a;
    const T bb = b*b;
    const T cc = c*c;
    const T caca = cosalpha * cosalpha;
    const T cbcb = cosbeta * cosbeta;
    const T cgcg = cosgamma * cosgamma;
    const T O_O = abc * abc * ( T(1) - caca - cbcb - cgcg + ccc + ccc );
   
    matrix<T> ans( 3, 3 );
    ans[0][0] = bb * cc * ( T(1) - caca );
    ans[0][1] = abc * c * f(cosalpha, cosbeta, cosgamma);
    ans[0][2] = abc * b * f(cosgamma, cosalpha, cosbeta);
    ans[1][0] = ans[0][1];
    ans[1][1] = aa * cc * ( T(1) - cbcb );
    ans[1][2] = abc * a * f(cosbeta, cosgamma, cosalpha);
    ans[2][0] = ans[0][2];
    ans[2][1] = ans[1][2];
    ans[2][2] = aa * bb * ( T(1) - cgcg );
    ans /= O_O;

    return ans;
}

}//namespace feng

#endif//_RECIPROCAL_METRIC_TENSOR_HPP_INCLUDED_DFOSIJEOIUSDFLKJALDFOI4ELKJDLKJFSAINTOSOIJTLKASOFIDJ3LKIJHSFALKJFDL

