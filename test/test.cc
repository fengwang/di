#include <construct_a.hpp>

#include <iostream>
#include <algorithm>
#include <iterator>

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
    auto I = feng::construct_a<double>().make_i(100); 

    //feng::for_each( I.begin(), I.end(), [](double&d){ if (std::abs(d)<=1.0e-8) d=0; } );

    std::cout << "\nI=\n" <<  I;


    std::cout << "\ndiagonal value of I\n";
    std::copy( I.diag_begin(), I.diag_end(), std::ostream_iterator<double>(std::cout, "\t"));

    std::cout << "\nmiddle column of I\n";
    auto const mid = I.row() >> 1;
    std::copy( I.col_begin(mid), I.col_end(mid), std::ostream_iterator<double>(std::cout, "\n") );

    return 0;
}
