#ifndef _VOLUME_HPP_INCLUDED_FOSDIJ3098USFLKJ4098USFDLIJ49O8ILFDSIKJLVSFDOI3498SFDLIJLSFKDJALKJLSKJFCKMVNJKLDFLKJF
#define _VOLUME_HPP_INCLUDED_FOSDIJ3098USFLKJ4098USFDLIJ49O8ILFDSIKJLVSFDOI3498SFDLIJLSFKDJALKJLSKJFCKMVNJKLDFLKJF

#include <tri_ary.hpp>

#include <cmath>
#include <cassert>

namespace feng
{
    template<typename T>
    T volume( const tri_ary<T>& a, const tri_ary<T>& b, const tri_ary<T>& c )
    {
        return dot_product( a, crossproduct( b, c ) );
    }
    
    //Input:
    //      lattice parameters {a,b,c,\alpha,\beta,\gamma}
    //Output:
    //      volume of the lattice
    template<typename T>
    T volume( const T a, const T b, const T c, const T alpha, const T beta, const T gamma )
    {
        const T abc = a * b * c;
        const T ca = std::cos(alpha);
        const T cb = std::cos(beta);
        const T cg = std::cos(gamma);
        const T cacbcg = ca * cb * cg;
        const T det = abc * abc * ( T(1) - ca*ca - cb*cb - cg*cg + cacbcg + cacbcg );
        assert( det > T() );
        return std::sqrt( det );
    }

}//namespace feng

#endif//_VOLUME_HPP_INCLUDED_FOSDIJ3098USFLKJ4098USFDLIJ49O8ILFDSIKJLVSFDOI3498SFDLIJLSFKDJALKJLSKJFCKMVNJKLDFLKJF

