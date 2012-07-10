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

    auto A_(A);
    std::fill( A_.diag_begin(), A_.diag_end(), double(0));
    std::cout << "\nmax non-diag element of A real is ";
    double mx=0;
    feng::for_each( A_.begin(), A_.end(), [&mx](std::complex<double> const & c){ if (std::abs(c.real()) > mx) mx = std::abs(c.real()); } ); 
    std::cout << mx;

    std::cout << "\naccumulation of non-diag value A is ";
    double acc = 0;
    feng::for_each( A_.begin(), A_.end(), [&acc](std::complex<double> const & c){ acc += std::abs(c.real());} ); 
    std::cout << acc;

    //std::cout << "\nA=\n" << A << "\n";

    auto const offset = feng::construct_a<double>().make_gxy2_offset();

    feng::construct_a<double> :: matrix_type Is( A.row(), offset.col() );


    //std::cout << "\nThe gxy2_offset matrix is \n" << offset;

    unsigned long t1 = std::clock();

    for ( size_t i = 0; i != offset.col(); ++i )
    {
        //auto I = feng::construct_a<double>().make_new_i_with_offset(A, 100, offset.col_begin(i));
        auto I = feng::construct_a<double>(). make_new_i_with_offset(A,7.879, offset.col_begin(i));
        feng::for_each( I.begin(), I.end(), [](double&d){ if (std::abs(d)<=1.0e-8) d=0; } );
        auto const mid = I.row() >> 1;
        std::copy( I.col_begin(mid), I.col_end(mid), Is.col_begin(i) );

        std::cout << "\naccumulation of " << i << "th column is " << std::accumulate( I.col_begin(mid), I.col_end(mid), double(0));
    }

    unsigned long t2 = std::clock();

    std::cout << "\ntime consumed: " << t2 - t1 << "\n";

    //std::cout << Is;

    //std::cout << "\nTest Is is \n" << feng::construct_a<double>().make_is( 7.879 ); 




    return 0;
}

