#ifndef _LATTICE_HPP_INCLUDED_SFODIJSFLDKJ4O8IUSFLKDJASFOIDJSFD
#define _LATTICE_HPP_INCLUDED_SFODIJSFLDKJ4O8IUSFLKDJASFOIDJSFD

#include <matrix.hpp>

#include <cmath>

namespace feng
{

    template<typename T>
    struct lattice
    {
        typedef T value_type;

        value_type a;
        value_type b;
        value_type c;
        value_type alpha;
        value_type beta;
        value_type gamma;
    };//struct lattice


    template<typename T>
    const matrix<T> make_matrix( const lattice<T>& l )
    {
        const T a = l.a;
        const T b = l.b;
        const T c = l.c;
        const T ca = std::cos(l.alpha);
        const T cb = std::cos(l.beta);
        const T cg = std::cos(l.gamma);
        const T sg = std::sin(l.gamma);
        matrix<T> A(3, 3);
        A[0][0] = a;    A[0][1] = 0;                A[0][2] = 0;
        A[1][0] = b*cg; A[1][1] = b*sg;             A[1][2] = 0;
        A[2][0] = c*cb; A[2][1] = c*(ca-cb*cg)/sg;  A[2][2] = c * std::sqrt(T(1)-ca*ca-cb*cb-cg*cg+T(2)*ca*cb*cg) / sg;
        return A;
    }

}//namespace feng

#endif//_LATTICE_HPP_INCLUDED_SFODIJSFLDKJ4O8IUSFLKDJASFOIDJSFD

