#include <iostream>

#include <misc.hpp>
#include <crystal.hpp>

void f(){}

int main()
{
    double a = 0.5;
    double c = 1.0;
    const double pi = 3.14159265358979323846;
    double alpha = pi / 2.0;

    std::cout << "a tetragonal crystal with lattice parameters a=0.5, and c=1.0\n";
    std::cout << "It's direct metric tensor:\n";
    std::cout << feng::direct_metric_tensor(a,a,c,alpha,alpha,alpha) << std::endl;
    std::cout << "It's reciprocal metric tensor:\n";
    std::cout << feng::reciprocal_metric_tensor(a,a,c,alpha,alpha,alpha) << std::endl;

    //feng::matrix<double> p(1,3);
    //feng::matrix<double> q(1,3);
    feng::tri_ary<double> p( 1, 2, 0 );
    feng::tri_ary<double> q( 3, 1, 1 );

    auto rmt = feng::reciprocal_metric_tensor(a,a,c,alpha,alpha,alpha);

    auto theta = feng::included_angle( p, q, rmt );
    std::cout << "included angle of " << p << " and " << q << " is " << theta << "\n";

    feng::tri_ary<double> t1(1,1,0);
    feng::tri_ary<double> t2(1,1,1);

    std::cout << t1 << " X " << t2 << feng::cross_product( t1, t2 )  << "\n";

    std::cout << "\ndirect structure matrix is:" << feng::direct_structure_matrix(a,a,c,alpha,alpha,alpha);
    std::cout << "\nreciprocal structure matrix is:" << feng::reciprocal_structure_matrix(a,a,c,alpha,alpha,alpha);

    feng::crystal<double> C(a,a,c,alpha,alpha,alpha);
    feng::tri_ary<double> t3(2,3,1);
    std::cout << "\ndirect components of vector " << t3 << " is " << C.direct_to_cartesian(t3) << "\n";
    feng::tri_ary<double> t4(2,2,1);
    std::cout << "\ncartesian components of vector " << t4 << " is " << C.cartesian_to_reciprocal(t4) << "\n";
    std::cout << "\nreciprocal components of vector " << t4 << " is " << C.reciprocal_to_cartesian(t4) << "\n";

    feng::crystal<long double> C_Al( 0.4049, 0.4049, 0.4049, alpha, alpha, alpha );
    feng::tri_ary<long double> g220(2,2,0);
    std::cout << "\nreciprocal space length of " << g220 << " is " << C_Al.reciprocal_space_length(g220);

    feng::crystal<long double> C_Si( 5.43, 5.43, 5.43, pi/2.0, pi/2.0, pi/2.0 );
    std::cout << "\nThe direct metric tensor matrix of Si is \n" << C_Si.direct_metric_tensor() << "\n";
    std::cout << "\nThe reciprocal metric tensor matrix of Si is \n" << C_Si.reciprocal_metric_tensor() << "\n";
    std::cout << "\nThe reciprocal structure matrix of Si is \n" << C_Si.reciprocal_structure_matrix() << "\n";





    return 0;
}

