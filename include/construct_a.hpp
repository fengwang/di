#ifndef _CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF
#define _CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

#include <array.hpp>
#include <lattice.hpp>
#include <matrix.hpp>

#include <vector>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>

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
            auto const gmax = make_gmax();
            auto const R    = make_reverse_matrix();
            auto const Rg   = R * gx;
            auto const nr   = Rg.norm();
            return std::ceil(gmax/nr);

        }

        size_type make_beamy_max() const 
        {
            auto const gmax = make_gmax();
            auto const R    = make_reverse_matrix();
            auto const gy   = make_gyvec();
            auto const Rg   = R * gy;
            auto const nr   = Rg.norm();
            return std::ceil(gmax/nr);
        }

        const matrix_type make_tbeam() const 
        {
            auto const Bx = make_beamx_max();
            auto const By = make_beamy_max();
            auto const gy = make_gyvec();
            
            make_gyvec tbeams{Bx*By+Bx+By, 5};
            size_type index = 0;
            for(size_type i = 0; i != Bx+1; ++i )
                for(size_type j = 0; j != By+1; ++j )
                {
                   if ( 0 == i && 0 == j ) continue;
                   tbeams[index][3] = i;
                   tbeams[index][4] = j;
                   auto const vec = -i*gx -j*gy;
                   std::copy( vec.begin(), vec.end(), tbeams.row_begin(index++) );
                }

            return tbeams;
        }

        const matrix_type gaussian_electron( const matrix_type&/*z*/,
                                             const matrix_type& s,
                                             const value_type v
                                           ) const 
        {
            const value_type* const A = {0.3626, 0.9737, 2.7209, 1.7660};
            const value_type* const B = {0.4281, 3.5770, 19.3905, 64.3334};
            const value_type        D = 0.0044;
            auto ans = s;
            auto const factor = v/value_type(511) + value_type(1);
            std::for_each( ans.begin(), ans.end(), [A,B,factor] (value_type& s_)
                           {
                                value_type tmp[4];
                                std::fill(tmp, tmp+4, value_type());
                                auto const ss = s_ * s_;
                                std::transform(B, B+4, tmp, [ss](value_type b) { return std::exp(-b*ss); });
                                s_ = factor * std::inner_product(A, A+4, tmp, value_type());
                           }
                          );
            return ans;
        }

        const matrix_type calc_ug(const matrix_type& g, //reciprocal lattice vector matrix
                                  const matrix_type& M, //
                                  const matrix_type& Z, //Z number matrix
                                  const matrix_type& A, //Atomic position matrix
                                  const matrix_type& D, //list of corresponding DW factor
                                  const value_type   v  //acceleration voltage
                                 )const
        {
            assert( 3 == M.row() );
            assert( 3 == M.col() );
            auto const R = M.inverse();
            auto const S = g*R;             
            auto const r = S.row();
            std::vector<value_type> gr( r );
            matrix_type s( 1, r );
            for ( size_type i = 0; i < r; ++ i )
                s[0][i] = value_type(0.5) * std::sqrt(std::inner_product(S.row_begin(i), S.row_end(i), S.row_begin(i), value_type(0)));
            auto const omega = feng::inner_product( array_type(M[0][0], M[1][0], M[2][0]),
                               feng::cross_product( array_type(M[0][1], M[1][1], M[2][1]), 
                                                    array_type(M[0][2], M[1][2], M[2][2])) );
            auto const atomcellfacte = gaussian_electron(Z, s, v);
        }




        





























    };//sturct construct_a

}//namespace feng

#endif//_CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

