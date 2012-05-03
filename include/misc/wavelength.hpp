#ifndef _WAVELENGTH_HPP_INCLUDED_DOSIFJ4EOIINTOTHYEFKSIWORKISIFDJBETTERSETTINGSISWRONGPROBABILIYTIJHSFLKJDFLKJFKJFD
#define _WAVELENGTH_HPP_INCLUDED_DOSIFJ4EOIINTOTHYEFKSIWORKISIFDJBETTERSETTINGSISWRONGPROBABILIYTIJHSFLKJDFLKJFKJFD

#include <cmath>

namespace feng
{
    template<typename T>
    T wavelength ( const T acceleration_voltage )
    {
        static const T h = 6.6260695729e-34; //Plank constant
        static const T m = 9.1093829140e-15; //electron mass
        static const T e = 1.60217656535e-19;//elementary charge
        static const T c = 299792458;        //speed of light in vacuum
        const T E = acceleration_voltage;
        const T eE = e * E;
        const T meE = m * eE;
        const T mcc = m * c * c;
        return h * std::sqrt ( ( meE + meE ) * ( T ( 1 ) + eE / ( mcc + mcc ) ) );
    }

}//namespace feng

#endif//_WAVELENGTH_HPP_INCLUDED_DOSIFJ4EOIINTOTHYEFKSIWORKISIFDJBETTERSETTINGSISWRONGPROBABILIYTIJHSFLKJDFLKJFKJFD

