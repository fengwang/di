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
#include <complex>
#include <vector>
#include <iterator>

#include <iostream>

namespace feng
{
    template<typename T = double>
    struct construct_a
    {
        typedef std::size_t                     size_type;
        typedef T                               value_type; 
        typedef std::complex<value_type>        complex_type; 
        typedef array<value_type>               array_type;
        typedef matrix<value_type>              matrix_type;
        typedef matrix<complex_type>            complex_matrix_type;
        typedef construct_a                     self_type;
        typedef lattice<value_type>             lattice_type;

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

        construct_a()
        {
            kt      = array_type{0, 0, 0};
            zone    = array_type{1, 1, 0};
            gx      = array_type{0, 0, 1};
            v0      = 200.0;
            la      = lattice_type{5.43, 5.43, 5.43, 1.5707963268, 1.5707963268, 1.5707963268 };
        }

        const array_type make_gyvec() const 
        {
            return cross_product(zone, gx);
        }

        const matrix_type make_matrix() const 
        {
            const T a = la.a;
            const T b = la.b;
            const T c = la.c;
            const T ca = std::cos(la.alpha);
            const T cb = std::cos(la.beta);
            const T cg = std::cos(la.gamma);
            const T sg = std::sin(la.gamma);
            matrix<T> A(3, 3);
            A[0][0] = a;    A[0][1] = 0;                A[0][2] = 0;
            A[1][0] = b*cg; A[1][1] = b*sg;             A[1][2] = 0;
            A[2][0] = c*cb; A[2][1] = c*(ca-cb*cg)/sg;  A[2][2] = c * std::sqrt(T(1)-ca*ca-cb*cb-cg*cg+T(2)*ca*cb*cg) / sg;
            return A;
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

        const matrix_type make_tbeams() const 
        {
            auto const Bx = make_beamx_max();
            auto const By = make_beamy_max();
            auto const gy = make_gyvec();
            
            matrix_type tbeams{Bx*By+Bx+By, 5};
            size_type index = 0;
            for(size_type j = 0; j != By+1; ++j )
                for(size_type i = 0; i != Bx+1; ++i )
                {
                   if ( 0 == i && 0 == j ) continue;
                   tbeams[index][3] = i;
                   tbeams[index][4] = j;
                   auto const vec = -(i*gx) - (j*gy);
                   std::copy( vec.begin(), vec.end(), tbeams.row_begin(index++) );
                }

            return tbeams;
        }

        const matrix_type make_gaussian_electron( const matrix_type& s,
                                                  const value_type v
                                                ) const 
        {
            //const value_type        D = 0.0044;
            auto ans = s;
            auto const factor = v/value_type(511) + value_type(1);
            std::for_each( ans.begin(), ans.end(), [factor] (value_type& s_)
                           {
                                const value_type A[] = {0.3626, 0.9737, 2.7209, 1.7660};
                                const value_type B[] = {0.4281, 3.5770, 19.3905, 64.3334};
                                value_type tmp[4];
                                std::fill(tmp, tmp+4, value_type());
                                auto const ss = s_ * s_;
                                std::transform(B, B+4, tmp, [ss](value_type b) { return std::exp(-b*ss); });
                                s_ = factor * std::inner_product(A, A+4, tmp, value_type());
                           }
                          );
            ans = ans && ans;
            ans = ans && ans;
            ans = ans && ans;

            return ans;
        }

        const complex_matrix_type make_ug(const matrix_type& g, //reciprocal lattice vector matrix
                                          const matrix_type& M, //
                                          const matrix_type& A, //Atomic position matrix
                                          const matrix_type& D //list of corresponding DW factor
                                         )const
        {
            assert( 3 == M.row() );
            assert( 3 == M.col() );
            auto const R = M.inverse();
            matrix_type const G(g, range(0, g.row()), range(0, 3));
            auto const S = G*R;             
            auto const r = S.row();
            std::vector<value_type> gr( r );
            matrix_type s( 1, r );
            for ( size_type i = 0; i < r; ++ i )
                s[0][i] = value_type(0.5) * std::sqrt(std::inner_product(S.row_begin(i), S.row_end(i), S.row_begin(i), value_type(0)));
            auto const piomega =    3.141592553590 * feng::inner_product( array_type(M[0][0], M[1][0], M[2][0]),
                                                     feng::cross_product( array_type(M[0][1], M[1][1], M[2][1]), 
                                                                          array_type(M[0][2], M[1][2], M[2][2])) );
            auto const atomcellfacte = make_gaussian_electron(s, v0);

            auto ss = s;
            std::for_each(ss.begin(), ss.end(), [](value_type& v){v*=v;});
            complex_matrix_type dwss = D * ss;
            complex_matrix_type piag = A * G.transpose();
            auto fact = - dwss - piag * complex_type(0, 6.2831853071796);
            std::for_each(fact.begin(), fact.end(), [](complex_type& c){c=std::exp(c);});

            assert(fact.row() == atomcellfacte.row());
            assert(fact.col() == atomcellfacte.col());
            std::transform(fact.begin(), fact.end(), atomcellfacte.begin(), fact.begin(),
                           [piomega](const complex_type f,  const value_type a){return f*a/piomega;});

            complex_matrix_type Ug(fact.col(), 1);
            for ( size_type i = 0; i < fact.col(); ++i )
                Ug[i][0] = std::accumulate(fact.col_begin(i), fact.col_end(i), complex_type());

            return Ug;
        }

        const matrix_type make_dw_factor() const 
        {
            matrix_type m(8,1);
            std::fill(m.begin(), m.end(), 0.0044);
            return m;
        }

        const matrix_type make_atomic_positions() const 
        {
            const value_type pos_x[] = {-0.125, -0.125, 0.375, 0.375, 0.125, 0.125, 0.625, 0.625};
            const value_type pos_y[] = {-0.125, 0.375, -0.125, 0.375, 0.125, 0.625, 0.125, 0.625};
            const value_type pos_z[] = { -0.25,  0.25,   0.25, -0.25,   0.0,   0.5,   0.5, 0.0};
            matrix_type m(8,3);
            std::copy( pos_x, pos_x+8, m.col_begin(0) );
            std::copy( pos_y, pos_y+8, m.col_begin(1) );
            std::copy( pos_z, pos_z+8, m.col_begin(2) );
            return m;
        }

        struct complex_size
        {
            complex_type c;
            size_type    s;
            complex_size( const complex_type& c_ = complex_type(0,0), const size_type s_ = size_type(0) ) : 
                c(c_), s(s_) {}
        };//complex_size
        
        const matrix_type make_beams(const matrix_type& tbeams, const complex_matrix_type& Ug) const
        {
            assert( tbeams.row() == Ug.row() );
            assert( tbeams.col() == 5 );
            assert( Ug.col() == 1 );
            
            std::vector<complex_size> m( tbeams.row() );
            for ( size_type i = 0; i < m.size(); ++i )
                m[i] = complex_size(Ug[i][0], i);
            std::sort( m.begin(), m.end(), [](const complex_size& lhs, const complex_size& rhs){ return std::abs(lhs.c) >= std::abs(rhs.c);} );

            matrix_type beams{1, 3};
            auto const gy = make_gyvec();
            for ( size_type i = 0; i < m.size(); ++i )
            {
                const int ix = tbeams[m[i].s][3];
                const int iy = tbeams[m[i].s][4];
                matrix_type b{2,3};
                if ( ix )
                {
                    array_type const b0 = ix * gx + iy * gy;
                    array_type const b1 = - b0;
                    b[0][0] = b0[0]; b[0][1] = b0[1]; b[0][2] = b0[2];
                    b[1][0] = b1[0]; b[1][1] = b1[1]; b[1][2] = b1[2];
                    beams  = beams && b;
                }
                if ( iy )
                {
                    array_type const b0 = - ix * gx + iy * gy;
                    array_type const b1 = - b0;
                    b[0][0] = b0[0]; b[0][1] = b0[1]; b[0][2] = b0[2];
                    b[1][0] = b1[0]; b[1][1] = b1[1]; b[1][2] = b1[2];
                    beams  = beams && b;
                }
            }

            return beams;
        }
       
        const matrix_type  make_gd( const matrix_type& beams ) const 
        {
            matrix_type gd{beams.row(), beams.col()};
            auto const nm = beams.row() >> 1;
            for ( size_t i = 1; i <= nm; ++i )
            {
                std::copy( beams.row_begin(i+i-1), beams.row_end(i+i-1), gd.row_begin(nm+i) );
                std::copy( beams.row_begin(i+i), beams.row_end(i+i), gd.row_begin(nm-i) );
            }
            return gd;
        }
       
        void make_gm( const matrix_type& gd, matrix_type& Gm0, 
                      matrix_type& Gm1, matrix_type& Gm2) const 
        {
            assert( gd.col() == 3 );
            auto const N = gd.row();
            const matrix_type gd0{ gd, range(0, N), range(0, 1) };
            const matrix_type gd1{ gd, range(0, N), range(1, 2) };
            const matrix_type gd2{ gd, range(0, N), range(2, 3) };

            Gm0 = repmat(gd0, 1, N) - repmat(gd0.transpose(), N, 1);
            Gm1 = repmat(gd1, 1, N) - repmat(gd1.transpose(), N, 1);
            Gm2 = repmat(gd2, 1, N) - repmat(gd2.transpose(), N, 1);

            std::copy( gd.col_begin(0), gd.col_end(0), Gm0.diag_begin() );
            std::copy( gd.col_begin(1), gd.col_end(1), Gm1.diag_begin() );
            std::copy( gd.col_begin(2), gd.col_end(2), Gm2.diag_begin() );
        }
       
        const matrix_type make_unique_beams( const matrix_type& Gm0,
                                             const matrix_type& Gm1,
                                             const matrix_type& Gm2
                                           ) const 
        {
            auto const N = Gm0.row();

            array_type vec( Gm0[0][1], Gm1[0][1], Gm2[0][1] );
            auto const gy = make_gyvec();
            matrix_type bArray( 1, 7 );
            bArray[0][0] = 1;
            bArray[0][1] = 2;
            bArray[0][2] = inner_product(vec, gx);
            bArray[0][3] = inner_product(vec, gy);
            bArray[0][4] = Gm0[0][1];
            bArray[0][5] = Gm1[0][1];
            bArray[0][6] = Gm1[0][1];

            //one more row than matlab version, need checking
            matrix_type bArray_( 1, 7 );
            for ( size_type h = 0; h != N; ++h )
            {
                for ( size_type k = 0; k != N-h+1; ++k )
                {
                    if ( k == h ) continue;

                    vec[0] = Gm0[h][k];
                    vec[1] = Gm1[h][k];
                    vec[2] = Gm2[h][k];
                    
                    bool flag = false;
                    for ( size_type j = 0; j != bArray.row(); ++j  )
                        if ( vec == array_type{bArray[j][4], bArray[j][5], bArray[j][6]})
                        {
                            flag = true; break;
                        }

                    if (flag) continue;
                    bArray_[0][0] =h+1;
                    bArray_[0][1] =k+1;
                    bArray_[0][2] =inner_product(vec, gx);
                    bArray_[0][3] =inner_product(vec, gy);
                    bArray_[0][4] =Gm0[h][k];
                    bArray_[0][5] =Gm1[h][k];
                    bArray_[0][6] =Gm2[h][k];
                    bArray = bArray && bArray_;
                }//k loop
            }//h loop

            return bArray;
        }

        const complex_matrix_type make_ug( const matrix_type& bArray ) const 
        {
            auto const g = matrix_type{bArray, range{0, bArray.row()}, range{4,7}};
            auto const R = make_matrix();
            auto const A = make_atomic_positions();
            auto const D = make_dw_factor();
            return make_ug( g, R, A, D );
        }

        const complex_matrix_type operator ()() const 
        {
        
        }

    };//sturct construct_a

}//namespace feng

#endif//_CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

