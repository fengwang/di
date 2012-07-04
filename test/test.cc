//#include <construct_a.hpp>
#include <simulation.hpp>

#include <iostream>
#include <algorithm>
#include <iterator>

#include <ctime>

int main()
{
    //std::cout << feng::construct_a<double>()();
    //std::cout << feng::construct_a<double>().make_a();
    //std::cout << "\n";
    //std::cout << feng::construct_a<double>()(10);
    //std::cout << feng::construct_a<double>().make_i(10);
    //feng::construct_a<double>().make_eigen();

    //auto A = feng::construct_a<double>().make_a();

    //std::cout << feng::construct_a<double>().make_i_using_eigen_method(10);

    //std::cout << feng::expm(A);

    //auto I = feng::construct_a<double>().make_i(7.879); 
    //auto I = feng::construct_a<double>().make_i(100); 

    auto A = feng::construct_a<double>().make_new_a();
    feng::for_each( A.begin(), A.end(), [](std::complex<double>&d){ if(std::abs(d.real())<1.8e-8) d.real(0); if(std::abs(d.imag())<1.0e-8) d.imag(0); } );

    auto const offset = feng::construct_a<double>().make_gxy2_offset();

    feng::construct_a<double> :: matrix_type Is( A.row(), offset.col() );

    for ( size_t i = 0; i != offset.col(); ++i )
    {
        //auto I = feng::construct_a<double>().make_new_i(A,100);
        auto I = feng::construct_a<double>().make_new_i_with_offset(A,100, offset.col_begin(i));
        feng::for_each( I.begin(), I.end(), [](double&d){ if (std::abs(d)<=1.0e-8) d=0; } );
        auto const mid = I.row() >> 1;
        //std::cout << "\naccumulate middle column, we get:\n" << std::accumulate( I.col_begin(mid), I.col_end(mid), double(0) );
        //std::copy( I.col_begin(mid), I.col_end(mid), std::ostream_iterator<double>(std::cout, "\t"));
        //std::cout << "\n";

        std::copy( I.col_begin(mid), I.col_end(mid), Is.col_begin(i) );
    }

    std::cout << Is;

    return 0;
}
