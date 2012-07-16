#ifndef _XOVER_MANAGER_HPP_INCLUDED_SDFOIJROTIHASFO8498YASFIUHVBKJ9348T7YHSDKGJN34987UHSDGKJNBIUH489
#define _XOVER_MANAGER_HPP_INCLUDED_SDFOIJROTIHASFO8498YASFIUHVBKJ9348T7YHSDKGJN34987UHSDGKJNBIUH489

#include <matrix.hpp>
#include <gray_to_uint64.hpp>
#include <uint64_to_gray.hpp>


#include <cmath>
#include <vector>
#include <cstdint>


namespace ga
{
    //      returns an id that is selected
    // usage:
    //      auto& cs = feng::singleton<xover_manager>::instance();
    //      cs.initialize(1080); //should only in ga manager
    //      auto const selected_chromosome = cs(); 
    struct xover_manager
    {
        std::vector<double> weigh_array;

        xover_manager( const std::size_t n = 196 )
        {
            initialize( n );
        }

        void initialize( const std::size_t n = 196 )
        {
            weigh_array.resize(n);
            const double alpha = std::log(n) / (n+n-2);
            const double factor = std::exp(-alpha);
            double current = factor;

            for ( std::size_t i = 0; i != n; ++i )
            {
                weigh_array[i] = current;//exp(-(i+1)\alpha)
                current *= factor; 
            }

            // weigh array elements should look like this
            // { [a], [a+a^2], ..., [a+a^2+...+a^N] }
            for ( std::size_t i = 1; i != n; ++i )
            {
               weigh_array[i] += weigh_array[i-1]; 
            }

            auto const acc = *(weigh_array.rbegin());

            // weigh array elements should look like this
            // { [x], [xx], ..., [1] }
            std::for_each( weigh_array.begin(), weigh_array.end(), [acc](double& v){ v/=acc; } );
        }
            
        std::size_t operator()() const 
        {
            // a random number U[0,1]
            auto& vg = feng::singleton<vg::variate_generator<double>>::instance();
            auto const p = vg();

            return std::distance( weigh_array.begin(), std::upper_bound( weigh_array.begin(), weigh_array.end(), p ) ) - 1;
        }
    };//struct xover_manager 


    struct uniform_binary_xover
    {
        typedef std::uint64_t uint_type;
        void operator()( const uint_type father, const uint_type mother, uint_type& son, uint_type& daughter ) const 
        {
            uint_type const f = uint64_to_gray()(father);
            uint_type const m = uint64_to_gray()(mother);
            uint_type const lower_mask = 0x5555555555555555UL; 
            uint_type const upper_mask = 0xaaaaaaaaaaaaaaaaUL; 

            son      = gray_to_uint64()( (f & lower_mask) | (m & upper_mask) );
            daughter = gray_to_uint64()( (m & lower_mask) | (f & upper_mask) );
        }
    
    };//struct binary_xover 

    // usage:
    // feng::matrix<uint64_t> mf;
    // feng::matrix<uint64_t> mm;
    // feng::matrix<uint64_t> ms;
    // feng::matrix<uint64_t> md;
    // //crossover mf and mm, result in ms and md
    // uniform_binary_matrix_xover()( mf, mm, ms, md );
    struct uniform_binary_matrix_xover
    {
        template<typename Symmetric_Binary_Matrix_Type>
        void
        operator()( const Symmetric_Binary_Matrix_Type& father, const Symmetric_Binary_Matrix_Type& mother, Symmetric_Binary_Matrix_Type& son, Symmetric_Binary_Matrix_Type& daughter ) const 
        {
            typedef typename Symmetric_Binary_Matrix_Type::value_type value_type;
            assert( father.row() == father.col() );
            assert( mother.row() == father.col() );
            assert( mother.row() == mother.col() );

            std::size_t n = mother.row();
            son.resize(n, n);
            daughter.resize(n, n);

            for ( std::size_t i = 1; i != n; ++i )
            {
                feng::for_each( father.upper_diag_begin(i), father.upper_diag_end(i), 
                                mother.upper_diag_begin(i),
                                son.upper_diag_begin(i),
                                daughter.upper_diag_begin(i),
                                []( const value_type f, const value_type m, value_type& s, value_type& d )
                                {
                                    uniform_binary_xover()( f, m, s, d );
                                } 
                             );
                //make symmetric
                std::copy( son.upper_diag_begin(i), son.upper_diag_end(i), son.lower_diag_begin(i) );
                std::copy( daughter.upper_diag_begin(i), daughter.upper_diag_end(i), daughter.lower_diag_begin(i) );
            }
            
        }
    };





};//namespace ga

#endif//_XOVER_MANAGER_HPP_INCLUDED_SDFOIJROTIHASFO8498YASFIUHVBKJ9348T7YHSDKGJN34987UHSDGKJNBIUH489
