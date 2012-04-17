#ifndef _GCD_HPP_INCLUDED_FSOIDJ4398UFAKLSJSFLKJNVXMCNSDFKJHNASLKFJOPIEURJLKAFJOIERLKGDSLFKSSKJFLFLKDJSFLDSKJ
#define _GCD_HPP_INCLUDED_FSOIDJ4398UFAKLSJSFLKJNVXMCNSDFKJHNASLKFJOPIEURJLKAFJOIERLKGDSLFKSSKJFLFLKDJSFLDSKJ

#include <cstddef>
#include <algorithm>

namespace feng
{

    template<typename T = unsigned long long>
    T gcd ( T u, T v )
    {
        std::size_t shift;

        if ( u == 0 || v == 0 )
        { return u | v; }

        for ( shift = 0; ( ( u | v ) & 1 ) == 0; ++shift )
        {
            u >>= 1;
            v >>= 1;
        }

        while ( ( u & 1 ) == 0 ) { u >>= 1; }

        do
        {
            while ( ( v & 1 ) == 0 ) { v >>= 1; }

            if ( u > v ) { std::swap ( u, v ); }

            v = v - u;
        }
        while ( v != 0 );

        return u << shift;
    }

}//namespace feng

#endif//_GCD_HPP_INCLUDED_FSOIDJ4398UFAKLSJSFLKJNVXMCNSDFKJHNASLKFJOPIEURJLKAFJOIERLKGDSLFKSSKJFLFLKDJSFLDSKJ
