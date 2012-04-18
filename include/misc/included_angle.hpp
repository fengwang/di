#ifndef _INCLUDED_ANGLE_HPP_INCLUDED_OFSIJO48IULFSDKJSDFLKJLSDFKAJSFADLKIJ9438IUELKASDFJASSFDIUJLVKCXJ8IUF98U98FUY
#define _INCLUDED_ANGLE_HPP_INCLUDED_OFSIJO48IULFSDKJSDFLKJLSDFKAJSFADLKIJ9438IUELKASDFJASSFDIUJLVKCXJ8IUF98U98FUY

#include <misc/dot_product.hpp>
#include <matrix.hpp>
#include <tri_ary.hpp>

#include <cassert>
#include <cmath>

namespace feng
{
#if 0
               -1( dot(p,q) ) 
    theta = cos  (----------)
                 ( |p| |q|  ) 
#endif

template<typename T, std::size_t D, typename A>
T
included_angle( const tri_ary<T>& p, const tri_ary<T>& q, const matrix<T,D,A>& g )
{
    assert( 3 == g.row() );
    assert( 3 == g.col() );

    const T pq = dot_product( p, q, g );
    const T pp = dot_product( p, p, g );
    const T qq = dot_product( q, q, g );

    return std::acos( std::sqrt(pq*pq/(pp*qq)) );
}

}//namespace feng

#endif//_INCLUDED_ANGLE_HPP_INCLUDED_OFSIJO48IULFSDKJSDFLKJLSDFKAJSFADLKIJ9438IUELKASDFJASSFDIUJLVKCXJ8IUF98U98FUY
