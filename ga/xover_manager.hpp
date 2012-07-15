#ifndef _XOVER_MANAGER_HPP_INCLUDED_SDFOIJROTIHASFO8498YASFIUHVBKJ9348T7YHSDKGJN34987UHSDGKJNBIUH489
#define _XOVER_MANAGER_HPP_INCLUDED_SDFOIJROTIHASFO8498YASFIUHVBKJ9348T7YHSDKGJN34987UHSDGKJNBIUH489

#include <cmath>
#include <vector>
#include <cstdint>

namespace ga
{
    // usage:
    //      auto& cs = feng::singleton<chromosome_selector>::instance();
    //      cs.initialize(1080); //should only in ga manager
    //      auto const selected_chromosome = cs(); 
    struct chromosome_selector
    {
        std::vector<double> weigh_array;

        chromosome_selector( const std::size_t n = 196 )
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
    };//struct chromosome_selector 





};//namespace ga

#endif//_XOVER_MANAGER_HPP_INCLUDED_SDFOIJROTIHASFO8498YASFIUHVBKJ9348T7YHSDKGJN34987UHSDGKJNBIUH489
