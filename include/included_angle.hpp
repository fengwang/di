#ifndef _INCLUDED_ANGLE_HPP_INCLUDED_OFSIJO48IULFSDKJSDFLKJLSDFKAJSFADLKIJ9438IUELKASDFJASSFDIUJLVKCXJ8IUF98U98FUY
#define _INCLUDED_ANGLE_HPP_INCLUDED_OFSIJO48IULFSDKJSDFLKJLSDFKAJSFADLKIJ9438IUELKASDFJASSFDIUJLVKCXJ8IUF98U98FUY

#include <matrix.hpp>

#include <dot_product.hpp>

#include <cassert>
#include <cmath>

namespace feng
{
#if 0
               -1( dot(p,q) ) 
    theta = cos  (----------)
                 ( |p| |q|  ) 
#endif
template<typename T>
T
included_angle( const matrix<T>& p, const matrix<T>& q, const matrix<T>& g )
{
    assert( 1 == p.row() ); 
    assert( 1 == q.row() ); 
    assert( 3 == p.col() );
    assert( 3 == q.col() );
    assert( 3 == g.row() );
    assert( 3 == g.col() );

    const T pq = dot_product( p, q, g );
    const T pp = dot_product( p, p, g );
    const T qq = dot_product( q, q, g );

    return std::acos( std::sqrt(pq*pq/(pp*qq)) );
}

}//namespace feng

#endif//_INCLUDED_ANGLE_HPP_INCLUDED_OFSIJO48IULFSDKJSDFLKJLSDFKAJSFADLKIJ9438IUELKASDFJASSFDIUJLVKCXJ8IUF98U98FUY
