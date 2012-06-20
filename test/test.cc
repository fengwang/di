#include <construct_a.hpp>

#include <iostream>


int main()
{
    //std::cout << feng::construct_a<double>()();
    //std::cout << feng::construct_a<double>().make_a();
    //std::cout << "\n";
    std::cout << feng::construct_a<double>()(10);
    //feng::construct_a<double>().make_eigen();

    //auto const A = feng::construct_a<double>().make_a();

    //std::cout << A << "\n";

    //std::cout << feng::construct_a<double>().make_i_using_eigen_method(10);

    //std::cout << feng::expm(A);

    return 0;
}
