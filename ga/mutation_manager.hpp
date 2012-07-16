#ifndef _MUATATION_MANAGER_HPP_INCLUDED_SDOI4E98UAFSKLJSAFOIHJASFKJHASFKHIUYHASFKJHYUHSFDKJHAIUHRFSD
#define _MUATATION_MANAGER_HPP_INCLUDED_SDOI4E98UAFSKLJSAFOIHJASFKJHASFKHIUYHASFKJHYUHSFDKJHAIUHRFSD

#include <vg.hpp>
#include <singleton.hpp>

#include <uint64_to_gray.hpp>
#include <gray_to_uint64.hpp>

#include <cstddef>
#include <cstdint>

namespace ga
{
    // usage:
    //          auto& mm = feng::singleton<mutation_manager>::instance();
    //          mm.initialize( 1234 );
    //          std::size_t the_mutation_chromosome_id = mm();
    struct mutation_manager
    {
        std::size_t n;

        mutation_manager( const std::size_t n_ = 196 ) : n(n_) {}

        void initialize( const std::size_t n_ )
        {
            n = n_;
        }

        std::size_t operator()() const 
        {
            auto& vg = feng::singleton<vg::variate_generator<double>>::instance();
            return vg() * n;
        }

        std::size_t operator()( const std::size_t n_ ) const 
        {
            auto& vg = feng::singleton<vg::variate_generator<double>>::instance();
            return vg() * n_;
        }
    };//struct mutation_manager 


    struct binary_mutation
    {
        typedef std::uint64_t uint_type;
        // 1) select a bit randomly
        // 2) flip the bit
        void operator()( uint_type& u ) const 
        {
            auto& vg = feng::singleton<vg::variate_generator<double>>::instance();
            uint_type const total_pos = sizeof(uint_type) << 3; 
            // 1)
            uint_type const mask_pos = total_pos * vg();
            uint_type const mask = 1 << mask_pos;
            u = uint64_to_gray()(u);
            // 2)
            u ^= mask;
            u = gray_to_uint64()(u);
        }
    };//struct binary_mutation 

    // usage:
    //      feng::matrix<uint64_t> m;
    //      binary_matrix_mutation()(m);
    struct binary_matrix_mutation
    {
        template<typename Symmetric_Binary_Matrix>
        void operator()( Symmetric_Binary_Matrix& sbm ) const 
        {
            assert( sbm.row() == sbm.col() );

            std::size_t const n = sbm.row();     
            for ( std::size_t i = 1; i != n; ++i )
            {
                std::transform( sbm.upper_row_begin(i), sbm.upper_row_end(i), sbm.upper_row_begin(i), binary_mutation() );
                std::copy( sbm.upper_row_begin(i), sbm.upper_row_end(i), sbm.lower_row_begin(i) );
            }
        }
    
    };



}//namespace ga

#endif//_MUATATION_MANAGER_HPP_INCLUDED_SDOI4E98UAFSKLJSAFOIHJASFKJHASFKHIUYHASFKJHYUHSFDKJHAIUHRFSD

