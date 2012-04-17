#ifndef _DOT_PRODUCT_HPP_INCLUDED_OFSIJD98U43TLKAFSJDLSAKFDJSLDKFJFSDLKJFSDALKJSFLKDJSLFDKJVSFOIJSAFOIJSALKIJIOJOII
#define _DOT_PRODUCT_HPP_INCLUDED_OFSIJD98U43TLKAFSJDLSAKFDJSLDKFJFSDLKJFSDALKJSFLKDJSLFDKJVSFOIJSAFOIJSALKIJIOJOII

#include <matrix.hpp>
#include <tri_ary.hpp>

#include <cassert>
#include <algorithm>
#include <numeric>

namespace feng
{

#if 0

   dot(p,q) = p_i q_{i,j} q_j

#endif
template<typename T, std::size_t D, typename A>
T
dot_product( const tri_ary<T>& p, const tri_ary<T>& q, const matrix<T,D,A>& g )
{
   assert(3==g.row());
   assert(3==g.col());

   return ( p * g ) * q;
}

}//feng

#endif//_DOT_PRODUCT_HPP_INCLUDED_OFSIJD98U43TLKAFSJDLSAKFDJSLDKFJFSDLKJFSDALKJSFLKDJSLFDKJVSFOIJSAFOIJSALKIJIOJOII

