#ifndef _PROBABILITY_MANAGER_HPP_INCLUDED_SFOJIWI8HAFKHJASFODIUH34987YAFOISUH3487YFEI8HVDKJVJKBASIUW
#define _PROBABILITY_MANAGER_HPP_INCLUDED_SFOJIWI8HAFKHJASFODIUH34987YAFOISUH3487YFEI8HVDKJVJKBASIUW

#include <time_manager.hpp>
#include <singleton.hpp>
#include <vg.hpp>

#include <cmath>

namespace ga
{
    //NOTE: Only two probability managers implemented here, so private inheritance method not employed

    // usage:
    //      auto& xpm = feng::singleton<xover_probability_manager>::instance();
    //      auto p = xpm(); //get the current probability
    template<typename T = double>
    struct xover_probability_manager
    {
        typedef T value_type;
        value_type alpha;

        //crossover rate should be very high at first, then as GA running, it drops to a low rate
        //here set the probability to 1 at the start time, 0.25 at the end time
        xover_probability_manager( const value_type alpha_ = 1.38629436111989061883 ) : alpha( alpha_ )
        {}

        value_type operator()() const
        {
            auto& tm = feng::singleton<time_manager>::instance();
            value_type elapse_t = tm.elapse_time();
            value_type planning_t = tm.planning_time();
            return std::exp( -alpha * elapse_t / planning_t );
        }

    };//struct xover_probability_manager

    // usage:
    //      auto& mpm = feng::singleton<mutation_probability_manager>::instance();
    //      auto p = mpm(); //current mutation probability
    template<typename T = double>
    struct mutation_probability_manager
    {
        typedef T value_type;
        value_type alpha;
        value_type beta;
        //mutation rate should be very high at first, then as GA running, it drops to a low rate
        //here set the probability to 0.25 at the start time, 0.01 at the end time
        mutation_probability_manager( value_type const alpha_ = 3.21887582486820074921,
                                      value_type const beta_ = 1.38629436111989061883 ) :
            alpha( alpha_ ), beta( beta_ )
        {}

        value_type operator()() const
        {
            auto& tm = feng::singleton<time_manager>::instance();
            value_type elapse_t = tm.elapse_time();
            value_type planning_t = tm.planning_time();
            return std::exp( -alpha * elapse_t / planning_t - beta );
        }

    };//struct mutation_probability_manager

}//namespace ga

#endif//_PROBABILITY_MANAGER_HPP_INCLUDED_SFOJIWI8HAFKHJASFODIUH34987YAFOISUH3487YFEI8HVDKJVJKBASIUW

