#ifndef _INTERPLANAR_SPACING_HPP_INCLUDED_FSODIJHSAFLDKJ4389IUSFLKJD90843LKFSSFUI4398KLIFSDJKLHVKFSJHO9REWI9483433
#define _INTERPLANAR_SPACING_HPP_INCLUDED_FSODIJHSAFLDKJ4389IUSFLKJD90843LKFSSFUI4398KLIFSDJKLHVKFSJHO9REWI9483433

#include <matrix.hpp>
#include <dot_product.hpp>

#include <cmath>

namespace feng
{

    #if 0
            Input:
                    g       reciprocal lattice vector
                    g_      reciprocal metric tensor
            Output:
                    the interplanar spacing d_{hkl}
    #endif

    template<typename T>
    T
    interplanar_spacing( const matrix<T>& g, const matrix<T>& g_ )
    {
        const T gg = dot_product( g, g, g_ );
        return T(1) / std::sqrt(gg);
    }

}//namespace feng


#endif//_INTERPLANAR_SPACING_HPP_INCLUDED_FSODIJHSAFLDKJ4389IUSFLKJD90843LKFSSFUI4398KLIFSDJKLHVKFSJHO9REWI9483433

