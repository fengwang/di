#include <construct_a.hpp>

#include <iostream>


int main()
{
    using namespace feng;
    using namespace std;

    construct_a<double> ca;

//    cout << ca.make_gyvec() << endl;
//    cout << ca.make_matrix() << endl;
//    cout << ca.make_reverse_matrix() << endl;
//    cout << ca.make_matrix() * ca.make_reverse_matrix() << endl;
//    cout << ca.make_gscale() << endl;
//    cout << ca.make_gmax() << endl;
//    cout << ca.make_beamx_max() << endl;
//    cout << ca.make_beamy_max() << endl;
//    cout << ca.make_tbeams() << endl;
//    cout << ca.make_dw_factor() << endl;
//    cout << ca.make_atomic_positions() << endl;
    //auto const tbeams = ca.make_tbeams();
    //auto const Ug = ca.make_ug(tbeams, ca.make_matrix(), ca.make_atomic_positions(), ca.make_dw_factor());
    //auto const beams = ca.make_beams( tbeams, Ug );

    //auto const beams = ca.make_beams();
    auto const beams = ca.make_beam_vector();
    //auto const beams = ca.make_beam_array();
    auto const gd = ca.make_gd(beams);
    //matrix<double> gm0, gm1, gm2;
    //ca.make_gm( gd, gm0, gm1, gm2 );
    auto const gm = ca.make_gm(gd);
    //auto const ubeams = ca.make_unique_beams( gm0, gm1, gm2  );
    auto const ubeams = ca.make_unique_beams( gm  );
    //std::cout << ca.make_ug( ubeams );
    //std::cout << ubeams;

    return 0;
}
