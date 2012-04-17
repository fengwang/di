#include <iostream>

#include <direct_metric_tensor.hpp>
#include <dot_product.hpp>
#include <included_angle.hpp>
#include <delta.hpp>
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

/*
    p[0][0] = 1;
    p[0][1] = 2;
    p[0][2] = 0;
    q[0][0] = 3;
    q[0][1] = 1;
    q[0][2] = 1;
*/

    auto rmt = feng::reciprocal_metric_tensor(a,a,c,alpha,alpha,alpha);

    auto theta = feng::included_angle( p, q, rmt );
    std::cout << "included angle of " << p << " and " << q << " is " << theta << "\n";


    return 0;
}

