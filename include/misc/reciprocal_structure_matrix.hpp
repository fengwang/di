#ifndef _RECIPROCAL_STRUCT_MATRIX_HPP_INCLUDED_SOFIJOI4EULFASKJLKSFDJS8888888888888888888888888888888888SFDSFDSDFF
#define _RECIPROCAL_STRUCT_MATRIX_HPP_INCLUDED_SOFIJOI4EULFASKJLKSFDJS8888888888888888888888888888888888SFDSFDSDFF

#include <misc/direct_structure_matrix.hpp>

namespace feng
{
    template<typename T>
    const matrix<T>
    reciprocal_structure_matrix( const T a, const T b, const T c, const T alpha, const T beta, const T gamma )
    {
        return direct_structure_matrix( a, b,  c, alpha, beta, gamma ).inverse();
    }

}//namespace feng

#endif//_RECIPROCAL_STRUCT_MATRIX_HPP_INCLUDED_SOFIJOI4EULFASKJLKSFDJS8888888888888888888888888888888888SFDSFDSDFF

