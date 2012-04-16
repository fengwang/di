#ifndef _RECIPROCAL_VECTOR_HPP_INCLUDED_SFODIJ4O8IUFSALKJH4398UYAFSLKIJH4398IUYASFKLIJH98T4YLKAJFSHVKJMNSKAFJDHFSJK
#define _RECIPROCAL_VECTOR_HPP_INCLUDED_SFODIJ4O8IUFSALKJH4398UYAFSLKIJH4398IUYASFKLIJH98T4YLKAJFSHVKJMNSKAFJDHFSJK

#include <matrix.hpp>

#include <cassert>
#include <cmath>
#include <algorithm>

namespace feng
{

#if 0
        dot(a_i, a^*_j) = delta_{i,j}
#endif

template<typename T>
bool
reciprocal_vector( const matrix<T>& a, const matrix<T>& b, const matrix<T>& c,
                   matrix<T>& a_,      matrix<T>& b_,      matrix<T>& c_ )
{
    assert( 1 == a.row() );
    assert( 1 == b.row() );
    assert( 1 == c.row() );
    assert( 3 == a.col() );
    assert( 3 == b.col() );
    assert( 3 == c.col() );

    a_.resize(1,3);
    b_.resize(1,3);
    c_.resize(1,3);

    const matrix<T> M = a && 
                        b && 
                        c;
    const matrix<T> M_ = M.inverse();

    //check if the matrix is inversable
    for(auto m : M_ )
        if (isnan(m) || isinf(m))
            return false;

    std::copy(M_.col_begin(0), M_.col_end(0), a_.begin());
    std::copy(M_.col_begin(1), M_.col_end(1), b_.begin());
    std::copy(M_.col_begin(2), M_.col_end(2), c_.begin());

    return true;
}

}//namespace feng

#endif//_RECIPROCAL_VECTOR_HPP_INCLUDED_SFODIJ4O8IUFSALKJH4398UYAFSLKIJH4398IUYASFKLIJH98T4YLKAJFSHVKJMNSKAFJDHFSJK


