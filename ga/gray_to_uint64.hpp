#ifndef _GRAY_TO_INT64_HPP_INCLUDED_DSOFIJR98AFSIHJ98R3YAOWI3498YAOISFUHASKFDHJASIFUH3HU7AIUHFDKJVKJ
#define _GRAY_TO_INT64_HPP_INCLUDED_DSOFIJR98AFSIHJ98R3YAOWI3498YAOISFUHASKFDHJASIFUH3HU7AIUHFDKJVKJ

#include <cstdint>

namespace ga
{

    // usage:
    //          uint64_t n;
    //          auto i = gray_to_int_64()(n); //convert n to normal int
    // 
    //          matrix<uint64_t> m;
    //          feng::for_each( m.begin(), m.end(), [](uint64_t& v){ v = gray_to_uint64()(v); } );
    //          std::transform( m.begin(), m.end(), m.begin(), gray_to_uint64() );
    struct gray_to_uint64
    {
        typedef std::uint64_t value_type;

        const value_type operator()( const value_type v ) const
        {
            value_type num = v;
            unsigned int const numBits = 64;
            for ( unsigned int shift = 1; shift < numBits; shift *= 2 )
                num ^= num >> shift;
            return num;
        }

    };//gray_to_uint64
};//namespace ga

#endif//_GRAY_TO_INT64_HPP_INCLUDED_DSOFIJR98AFSIHJ98R3YAOWI3498YAOISFUHASKFDHJASIFUH3HU7AIUHFDKJVKJ

