#include <iostream>

#include <direct_metric_tensor.hpp>
#include <dot_product.hpp>
#include <included_angle.hpp>
#include <delta.hpp>
#include <miller.hpp>
#include <volume.hpp>
#include <reciprocal_vector.hpp>
#include <reciprocal_metric_tensor.hpp>
#include <interplanar_spacing.hpp>

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


    return 0;
}

