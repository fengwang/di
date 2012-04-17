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
reciprocal_vector( const tri_ary<T>& a, const tri_ary<T>& b, const tri_ary<T>& c,
                   tri_ary<T>& a_,      tri_ary<T>& b_,      tri_ary<T>& c_ )
{
    matrix<T> M(3,3);
    M[0][0] = a.x();
    M[0][1] = a.y();
    M[0][2] = a.z();
    M[1][0] = b.x();
    M[1][1] = b.y();
    M[1][2] = b.z();
    M[2][0] = c.x();
    M[2][1] = c.y();
    M[2][2] = c.z();
    const matrix<T> M_ = M.inverse();

    //check if the matrix is inversable
    for(auto m : M_ )
        if (isnan(m) || isinf(m))
            return false;

    a_.x() = M_[0][0];
    a_.y() = M_[1][0];
    a_.z() = M_[2][0];
    b_.x() = M_[0][1];
    b_.y() = M_[1][1];
    b_.z() = M_[2][1];
    c_.x() = M_[0][2];
    c_.y() = M_[1][2];
    c_.z() = M_[2][2];

    return true;
}

}//namespace feng

#endif//_RECIPROCAL_VECTOR_HPP_INCLUDED_SFODIJ4O8IUFSALKJH4398UYAFSLKIJH4398IUYASFKLIJH98T4YLKAJFSHVKJMNSKAFJDHFSJK


