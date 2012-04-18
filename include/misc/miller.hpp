#ifndef _MILLER_HPP_INCLUDED_FSDOIJ43908UASFLJLSFDDDDDDDDDDDDDDDDDDDDDDDDIEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#define _MILLER_HPP_INCLUDED_FSDOIJ43908UASFLJLSFDDDDDDDDDDDDDDDDDDDDDDDDIEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#include <tri_ary.hpp>
#include <quadri_ary.hpp>
#include <numeric/gcd.hpp>

#include <cmath>

namespace feng
{

    template<typename T>
    const tri_ary<T>
    miller_bravais_to_miller( const quadri_ary<T>& q )
    {
        return tri_ary<T>( q.u()-q.t(), q.v()-q.t(), q.w() );
    }

    template<typename T>
    const quadri_ary<T>
    miller_to_miller_vravais( const tri_ary<T>& q )
    {
        unsigned long long u = static_cast<unsigned long long>(std::abs( q.x()+q.x()-q.y() ));
        unsigned long long v = static_cast<unsigned long long>(std::abs( q.y()+q.y()-q.x() ));
        unsigned long long w = static_cast<unsigned long long>(std::abs( -q.x()-q.y() ));
        unsigned long long t = static_cast<unsigned long long>(std::abs( q.z()+q.z()+q.z() ));

        const unsigned long long g = gcd(gcd(u,v), gcd(w,t));
        return quadri_ary<T>(u/g, v/g, w/g, t/g);
    }

}//namespace feng



#endif//_MILLER_HPP_INCLUDED_FSDOIJ43908UASFLJLSFDDDDDDDDDDDDDDDDDDDDDDDDIEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

