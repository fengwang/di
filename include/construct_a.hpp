#ifndef _CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF
#define _CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

#include <array.hpp>
#include <lattice.hpp>
#include <matrix.hpp>

#include <cstddef>
#include <cmath>

namespace feng
{
    template<typename T = double>
    struct construct_a
    {
        typedef std::size_t         size_type;
        typedef T                   value_type; 
        typedef array<value_type>   array_type;
        typedef matrix<value_type>  matrix_type;
        typedef construct_a         self_type;
        typedef lattice<value_type> lattice_type;

        array_type      kt;
        array_type      zone;
        array_type      gx;
        value_type      v0;
        lattice_type    la;

        construct_a(const self_type&) = default;
        construct_a(self_type&&) = default;
        self_type& operator = (const self_type&) = default;
        self_type& operator = (self_type&&) = default;
        
        construct_a( const array_type& kt_, const array_type& zone_,
                     const array_type& gx_, const value_type  v0_,
                     const lattice_type& la_): 
                     kt(kt_), zone(zone_), gx(gx_), v0(v0_), la(la_)
        {}

        const array_type make_gyvec() const 
        {
            return cross_product(zone, gx);
        }

        const matrix_type make_matrix() const 
        {
            returnn ::make_matrix(la);
        }

        const matrix_type make_reverse_matrix() const 
        {
            return make_matrix().inverse();
        }
        
        const array_type make_gscale() const 
        {
            return array_type{value_type(1)/la.a, value_type(1)/la.b, value_type(1)/la.c}; 
        }

        const value_type make_gmax() const 
        {
            return value_type(1.25);
        }

        void plot_structure() const;

        size_type make_beamx_max() const 
        {
            auto gmax = make_gmax();
            auto R    = make_reverse_matrix();
            auto Rg   = R * gx;
            auto nr   = Rg.norm();
            return std::ceil(gmax/nr);
        }

        size_type make_beamy_max() const 
        {
            auto gmax = make_gmax();
            auto R    = make_reverse_matrix();
            auto gy   = make_gyvec();
            auto Rg   = R * gy;
            auto nr   = Rg.norm();
            return std::ceil(gmax/nr);
        }

        const matrix_type make_tbeam() const 
        {
            auto Bx = make_beamx_max();
            auto By = make_beamy_max();
            auto gy = make_gyvec();
            
            make_gyvec tbeams{Bx*By+Bx+By, 5};
            size_type index = 0;
            for(size_type i = 0; i != Bx+1; ++i )
                for(size_type j = 0; j != By+1; ++j )
                {
                   if ( 0 == i && 0 == j ) continue;
                   tbeams[index][3] = i;
                   tbeams[index][4] = j;
                   auto vec = -i*gx -j*gy;
                   std::copy( vec.begin(), vec.end(), tbeams.row_begin(index++) );
                }

            return tbeams;
        }





        





























    };//sturct construct_a

}//namespace feng

#endif//_CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

