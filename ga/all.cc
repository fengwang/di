#include<binary_matrix_to_real_complex_matrix.hpp>
#include<chromosome.hpp>
#include<gray_to_uint64.hpp>
#include<uint64_to_gray.hpp>
#include<population.hpp>
#include<probability_manager.hpp>
#include<singleton.hpp>
#include<time_manager.hpp>
#include<xover_manager.hpp>
#include<mutation_manager.hpp>
#include<random_matrix_initialize.hpp>

#include <genetic_algorithm.hpp>

#include <iostream>

int main()
{
    using namespace ga;
    //chromosome<feng::matrix<std::uint64_t>, symmetric_random_matrix_initialize, double> c(5,5);

    //std::cout << *(c.ch);

    //xover_selection_manager ss(256);

    //std::copy( ss.weigh_array.begin(), ss.weigh_array.end(), std::ostream_iterator<double>(std::cout, "\n"));



    //genetic_algorithm<> ga(100, 196);
    genetic_algorithm<> ga(200, 196);

    std::cout << "\n" << ga(25,25);

    //chromosome<> ch(9,9);

    //std::cout <<  *(ch.chrom);
    //std::cout <<  binary_matrix_to_real_complex_matrix_symmetric()(*(ch.chrom), 0.0, 0.06);




    return 0;
}

