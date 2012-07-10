#ifndef _GA_HPP_INCLUDED_DFSOIJ4983UASFDOIJASFKJAHSAFKDJDFGKJHSAFIUHJ498UHYASKFHGAISUFHFISUDHIUGHRKF
#define _GA_HPP_INCLUDED_DFSOIJ4983UASFDOIJASFKJAHSAFKDJDFGKJHSAFIUHJ498UHYASKFHGAISUFHFISUDHIUGHRKF

#include <vg.hpp>
#include <simulation.hpp>

#include <matrix.hpp>

#include <chrono>
#include <cstddef>
#include <cassert>
#include <complex>
#include <vector>
#include <memory>
#include <utility>
#include <cmath>
#include <cassert>

namespace feng
{

template<typename T=double>
struct chromosome
{
    typedef T value_type;
    typedef std::int32_t gene_value_type;
    typedef vg::variate_generator<value_type> vg_type;

    // a binary matrix representing the chromosome
    matrix<gene_value_type> gene_matrix;
    
    value_type residual;//|| I - I' ||^2
    value_type ratio; //probability to be selected
    bool is_a_fresh_man;//false if this residual for this chromosome has been calculated since last change in gene_matrix

    chromosome( const std::size_t n ) : gene_matrix(n, n)
    {
         //is there a c++11 singleton lib that accepts initializer arguments?
         vg_type& vg = feng::singleton<vg>::instance();
         for ( std::size_t i = 1; i != n; ++i )
         {
             //random gene
            feng::for_each( gene_matrix.upper_diag_begin(i), gene_matrix.upper_diag_end(i), [&](gene_value_type& g)
                            { g = static_cast<gene_value_type>( vg() * std::numeric_limits<gene_value_type>::max() ); } );
            //symmetric matrix
            std::copy( gene_matrix.upper_diag_begin(i), gene_matrix.upper_diag_end(i), gene_matrix.lower_diag_begin(i) );
         }
    }


};//chromosome

struct chromosome_to_a
{
    // translate chromosome to potential matrix A
    template<typename Chromosome_Type, typename Gxy2_Itor>
    matrix<std::complex<typename Chromosome_Type::value_type> > const 
    operator()( const Chromosome_Type& ch,  // chromosome to be translated
                Gxy2_Itor gxy2_begin, // including beams offsets
                const typename Chromosome_Type::value_type lower_bound = -0.06,  // min value
                const typename Chromosome_Type::value_type upper_bound = 0.06) const 
    {
        typedef typename Chromosome_Type::value_type value_type;
        typedef typename Chromosome_Type::gene_value_type gene_value_type;
        value_type const g_min = static_cast<value_type>(std::numeric_limits<gene_value_type>::min());
        value_type const g_max = static_cast<value_type>(std::numeric_limits<gene_value_type>::max());
        value_type const g_diff = g_max - g_min;
        value_type const bound_diff = upper_bound - lower_bound;

        assert( ch.gene_matrix.row() == ch.gene_matrix.col() );
        std::size_t const n = ch.gene_matrix.row();
        matrix<std::complex<value_type> > ans( n, n );

        for ( std::size_t i = 1; i != n; ++i )
        {
            std::for_each( ans.upper_diag_begin(i), ans.upper_diag_end(i), ch.gene_matrix.upper_diag_begin(i),
                           [g_min, bound_diff, g_diff, lower_bound](std::complex<value_type>& c, gene_value_type const g )
                           { 
                               const value_type vg = static_cast<value_type>(g);
                               c.real( (vg-g_min) * bound_diff / g_diff + lower_bound );
                           } );
            //symmetric matrix
            std::copy( ans.upper_diag_begin(i), ans.upper_diag_end(i), ans.lower_diag_begin(i) );
        }

        //update diagonal value of matrix A
        feng::for_each( ans.diab_begin(), ans.diag_end(), gxy2_begin, []( std::complex<value_type>& c, const value_type v ) { c.real(-v); } );

        return ans;
    }
};//chromosome_to_a 


struct residual_evaluator
{
    // calculate the residuals 
    template<typename A_Matrix, typename Offset_Itor, typename Current_Itor, typename T>
    T const 
    operator()( A_Matrix const & A,     //the scattering matrix A working
                Offset_Itor off_begin,  //introduced beam tilt
                Current_Itor current_begin, //target current
                T const pi_lambda_t = 7.879 ) const 
    {
        typedef typename Current_Itor::value_type value_type;
        auto I = feng::construct_a<value_type>().make_new_i_with_offset( A, pi_lambda_t, off_begin );
        T residual = 0;
        std::size_t const mid = I.row() >> 1;
        //the central column
        feng::for_each( I.col_begin(mid), I.col_end(mid), current_begin, [&residual](value_type const v1, value_type const v2)
                        {  auto const diff = v1 - v2; residual += diff*diff; } );
        return residual;
    }

};


#if 0
struct chromosome_residual_evaluator
{
    template<typename Chromosome_Type, typename Diagonal_Beam_Itor, typename Scattering_Current_Itor>
    typename Chromosome_Type::value_type 
    operator()( Chromosome_Type const & ch, 
                Diagonal_Beam_Itor beam_first, Diagonal_Beam_Itor beam_last,  //this itor should refer to make_new_a
                Scattering_Current_Itor current_first, Scattering_Current_Itor current_last ) const 
    {
        auto A = chromosome_to_ug()( ch ); 
        assert( std::distance(beam_first, beam_last) == A.row() );
        assert( std::distance(current_first, current_last) == A.row() );
        for ( std::size_t i = 0; i != A.row(); ++i )
            A[i][i] = *beam_first++;
    }

};//struct chromosome_fitness_evaluator 
#endif

template<typename T>
struct ga_manager
{
    typedef std::chrono::time_point<std::chrono::system_clock> time_point_type;
    typedef std::size_t size_type;
    typedef construct_a<T> simulation_type;
    typedef typename simulation_type::value_type value_type;
    typedef typename simulation_type::matrix_type matrix_type;
    typedef typename simulation_type::vector_type vector_type;
    typedef typename simulation_type::complex_vector_type complex_vector_type;
    typedef typename simulation_type::complex_matrix_type complex_matrix_type;
    typedef chromosome<T> chromosome_type;
    typedef std::shared_ptr<chromosome_type> chromosome_ptr_type;
    typedef std::vector<chromosome_ptr_type> chromosome_pool_type;

    //working chromosome pools

    chromosome_pool_type survival_pool;     //contains individuals that suvival after selection
    chromosome_pool_type father_xover_pool; //contains individuals selected for xover
    chromosome_pool_type mother_xover_pool; //contains individuals selected for xover
    chromosome_pool_type new_generation_pool; //contains individuals newly generated for xover pools
    chromosome_pool_type mutation_pool;     //contains individuals to be mutated
    chromosome_pool_type living_pool;       //contains individuals living now

    complex_matrix_type A; //make_a(); -- this matrix shall never be used, but for comparison with GA result
    matrix_type gxy2_offset; //make_gxy2_offset()
    vector_type gxy2;//make_gxy2()
    matrix_type Is; //make_is(pi_lambda_t)

    value_type pi_lambda_t;

    time_point_type start_time;
    size_type total_run_time;
    size_type elapse_time;

    ga_manager( value_type pi_lambda_t_ = 7.879,
                size_t const total_time = 1000) 
        :A(                 simulation_type().make_a()),
         gxy2_offset(       simulation_type().make_gxy2_offset()),
         gxy2(              simulation_type().make_gxy2()),
         Is(                simulation_type().make_is(pi_lambda_t_) ),
         pi_lambda_t(       pi_lambda_t_), 
         start_time(        std::chrono::system_clock::now()),
         total_run_time(    total_time), 
         elapse_time(       0)
    { }

    void update_elapse_time()
    {
        time_point_type const current_time = std::chrono::system_clock::now();
        elapse_time = std::chrono::duration_cast<std::chrono::seconds> (current_time - start_time).count();
    }

    bool is_time_out()
    {
        update_elapse_time();
        return elapse_time > total_run_time;
    }

    // crossover rate, from 1 to 0.25
    value_type xover_rate() 
    {
        update_elapse_time();
        value_type const alpha = 1.38;
        return std::exp( - alpha * elapse_time / total_run_time );
    }

    // mutation rate, from 0.5 to 0.001
    value_type mutation_rate() 
    {
        update_elapse_time();
        value_type const alpha = 0.69;
        return std::exp( - alpha * elapse_time / total_run_time ) - 0.5;
    }

    // individuals directly passed to next generation, from 0.25 to 0.75
    value_type survival_rate() 
    {
        update_elapse_time();
        value_type const alpha = -1.38;
        value_type const beta =  1.1;
        return std::exp( alpha + beta * elapse_time / total_run_time );
    }



};//struct ga_manager


}//namespace feng
#endif//_GA_HPP_INCLUDED_DFSOIJ4983UASFDOIJASFKJAHSAFKDJDFGKJHSAFIUHJ498UHYASKFHGAISUFHFISUDHIUGHRKF

