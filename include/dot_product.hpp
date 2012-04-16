#ifndef _DOT_PRODUCT_HPP_INCLUDED_OFSIJD98U43TLKAFSJDLSAKFDJSLDKFJFSDLKJFSDALKJSFLKDJSLFDKJVSFOIJSAFOIJSALKIJIOJOII
#define _DOT_PRODUCT_HPP_INCLUDED_OFSIJD98U43TLKAFSJDLSAKFDJSLDKFJFSDLKJFSDALKJSFLKDJSLFDKJVSFOIJSAFOIJSALKIJIOJOII

#include <matrix.hpp>

#include <cassert>
#include <algorithm>
#include <numeric>

namespace feng
{

#if 0

   dot(p,q) = p_i q_{i,j} q_j

#endif
template<typename T>
T
dot_product( const matrix<T>& p, const matrix<T>& q, const matrix<T>& g )
{
   assert(1==p.row()); 
   assert(1==q.row()); 
   assert(3==p.col());
   assert(3==q.col());
   assert(3==g.row());
   assert(3==g.col());

   const matrix<T> pg = p*g;
   return std::inner_product(pg.begin(), pg.end(), q.begin(), T());
}

}//feng

#endif//_DOT_PRODUCT_HPP_INCLUDED_OFSIJD98U43TLKAFSJDLSAKFDJSLDKFJFSDLKJFSDALKJSFLKDJSLFDKJVSFOIJSAFOIJSALKIJIOJOII

