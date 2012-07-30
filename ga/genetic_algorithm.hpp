#ifndef _GENETIC_ALGORITHM_HPP_INCLUDED_SFOJ3W94889YKLSFDJHXVCMSAJKL4389SFDJKSFAJKLSFAFIOOSFIDWREOFF
#define _GENETIC_ALGORITHM_HPP_INCLUDED_SFOJ3W94889YKLSFDJHXVCMSAJKL4389SFDJKSFAJKLSFAFIOOSFIDWREOFF

#include <memory>
#include <vector>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <thread>

#include <expm_evaluator.hpp>
#include <random_matrix_initialize.hpp>
#include <singleton.hpp>

namespace ga
{

    //chromosome<feng::matrix<std::uint64_t>, 
    //           symmetric_random_matrix_initialize, 
    //           double> 
    //          ch(13,13);
    template<
                class Chromosome_Representation_Type = feng::matrix<std::uint64_t>,
                class Chromosome_Representation_Initialization_Type = symmetric_random_matrix_initialize,
                class Fitness_Type = double,
                class Fitness_Comparision_Type=std::less<Fitness_Type>
            >
    struct chromosome
    {
        typedef typename Chromosome_Representation_Type::iterator iterator;
        typedef typename Chromosome_Representation_Type::const_iterator const_iterator;
        
        std::shared_ptr<Chromosome_Representation_Type>  chrom;
        Fitness_Type                                     fit;
        bool                                             modification_after_last_evaluation_flag;
        
        chromosome( const chromosome& other ): chrom(other.chrom), fit(other.fit), modification_after_last_evaluation_flag(other.modification_after_last_evaluation_flag)
        {}

        //creating is initialization
        template<typename ... Args>
        chromosome( Args ... args )
        {
            chrom = std::make_shared<Chromosome_Representation_Type>(args...);
            Chromosome_Representation_Initialization_Type()(*chrom);
            modification_after_last_evaluation_flag = true;
        }

        friend bool operator < ( const chromosome& lhs, const chromosome& rhs )
        {
            return Fitness_Comparision_Type()(lhs.fit, rhs.fit); 
        }

        friend bool operator == ( const chromosome& lhs, const chromosome& rhs )
        {
            return (lhs.chrom == rhs.chrom) || (lhs.fit == rhs.fit);
        }

        const_iterator begin() const
        {
            return (*chrom).begin();
        }

        const_iterator end() const
        {
            return (*chrom).end();
        }

        iterator begin()
        {
            return (*chrom).begin();
        }

        iterator end()
        {
            return (*chrom).end();
        }
    };//struct chromosome

    template<
                class Chromosome_Representation_Type = feng::matrix<std::uint64_t>,
                class Fitness_Type = double,
                class Problem_Representation_Type = feng::matrix<std::complex<Fitness_Type>>,
                class Chromosome_Problem_Translation_Method = binary_matrix_to_real_complex_matrix_symmetric,
                class Chromosome_Mutation_Method = binary_mutation,
                class Chromosome_Xover_Method = single_point_xover,
                class Evaluation_Method = feng::expm_evaluator<Fitness_Type>, 
                class Chromosome_Representation_Initialization_Type = symmetric_random_matrix_initialize,
                class Fitness_Comparision_Type = std::less<Fitness_Type>
            >
    struct genetic_algorithm
    {
        typedef chromosome<Chromosome_Representation_Type, Chromosome_Representation_Initialization_Type, Fitness_Type, Fitness_Comparision_Type> chromosome_type;
        typedef std::vector<chromosome_type> population_type; 

        std::size_t population_per_generation;
        
        population_type current_population;
        population_type selected_population_for_xover_father;
        population_type selected_population_for_xover_mother;

        Problem_Representation_Type pr;

        Fitness_Type best_fitness;
        Problem_Representation_Type best_one;

        // setup runtime and population for the genetic algorithm
        genetic_algorithm( const std::size_t runtime_in_second_=3600, const std::size_t population_per_generation_ = 256 ): population_per_generation(population_per_generation_)
        { 
            current_population.reserve( population_per_generation ); 
            selected_population_for_xover_father.reserve( population_per_generation ); 
            selected_population_for_xover_mother.reserve( population_per_generation ); 

            //Should all these managers move to private member list????
            //setup timer
            auto& t_manager = feng::singleton<time_manager>::instance();
            t_manager(runtime_in_second_);

            //setup selection manager
            auto& xs_manager = feng::singleton<xover_selection_manager>::instance();
            xs_manager.initialize( population_per_generation );

            //setup mutation manager
            auto& m_manager = feng::singleton<mutation_manager>::instance();
            m_manager.initialize( population_per_generation );
        }

        struct parallel_evaluator
        {
            template<typename Iterator, typename Representation>
            void operator()(Iterator first, Iterator last, Representation* pr_) const 
            {
                for ( auto pch =first; pch != last; ++pch )
                    if ( (*pch).modification_after_last_evaluation_flag )
                    {
                        Chromosome_Problem_Translation_Method()( *((*pch).chrom), *pr_ );
                        (*pch).fit = Evaluation_Method()(*pr_);
                        (*pch).modification_after_last_evaluation_flag = false;
                    }
            }
        };//struct parallel_evaluator 

        // ch_args... are just arguments for chromosome initialization
        template<typename ... Chromosome_Args>
        Problem_Representation_Type const 
        operator()( Chromosome_Args ... ch_args )  
        {
            //initialize first generation randomly
            for ( std::size_t i = 0; i != population_per_generation; ++i )
                current_population.push_back( chromosome_type( ch_args... ) );

            //std::cerr << "\ninitialized " << population_per_generation << " chromosomes.\n";


            std::size_t debug_counter = 0;
            best_fitness = std::numeric_limits<Fitness_Type>::max();
            srand ( unsigned ( time (NULL) ) ); //for random_shuffle
            Fitness_Type total_fit;

            Problem_Representation_Type pr_0;
            Problem_Representation_Type pr_1;
            Problem_Representation_Type pr_2;
            Problem_Representation_Type pr_3;

            std::ofstream elite_out( "elite.dat" );
            std::ofstream elite_fit( "elite_fit.dat");
            
            for (;;)
            {
            //evaluate
                //feng::for_each( current_population.begin(), current_population.end(), [](chromosome_type& chrom) { Evaluation_Method()(Chromosome_Problem_Translation_Method()(chrom)); });
                //std::cerr << "\nevaluating........";

#if 0
                for ( auto& ch : current_population )
                    if ( ch.modification_after_last_evaluation_flag )
                    {
                        Chromosome_Problem_Translation_Method()( *(ch.chrom), pr );
                        ch.fit = Evaluation_Method()(pr);
                        //Evaluation_Method()(Chromosome_Problem_Translation_Method()(ch));
                        ch.modification_after_last_evaluation_flag = false;
                    }
#endif
#if 1
                std::thread t0( parallel_evaluator(), current_population.begin(), current_population.begin()+(current_population.size()/4), &pr_0 );
                std::thread t1( parallel_evaluator(), current_population.begin()+(current_population.size()/4), current_population.begin()+(current_population.size()/2), &pr_1 );
                std::thread t2( parallel_evaluator(), current_population.begin()+(current_population.size()/2), current_population.begin()+(current_population.size()*3/4), &pr_2 );
                std::thread t3( parallel_evaluator(), current_population.begin()+(current_population.size()*3/4), current_population.end(), &pr_3 );
                //parallel_evaluator()( current_population.begin()+(current_population.size()*3/4), current_population.end(), &pr_3 );
                t0.join();
                t1.join();
                t2.join();
                t3.join();
#endif

                //using nth_element?
                //mutate duplicated chromosome?
                std::sort( current_population.begin(), current_population.end() );

                auto const elite_one = *(current_population.begin());
                Chromosome_Problem_Translation_Method()( *(elite_one.chrom), pr );

                //keep the best individual of the current generation
                if ( elite_one.fit < best_fitness )
                {
                    debug_counter++;

                    best_fitness = elite_one.fit;
                    best_one = pr;
                    
                    elite_out << best_one;
                    elite_fit << best_fitness;


                    if ( debug_counter & 0xf )
                    {
                        elite_out << "\n";
                        elite_fit << "\n";
                    }
                    else
                    {
                        elite_out << std::endl;
                        elite_fit << std::endl;
                    }
                }

                current_population.resize( std::distance( current_population.begin(), std::unique( current_population.begin(), current_population.end() )) );
                if ( current_population.size() < population_per_generation )
                {
                    //generate someone randomly and evaluate
                    for ( std::size_t i = current_population.size(); i != population_per_generation; ++i )
                    {
                        auto ch = chromosome_type( ch_args...);
                        Chromosome_Problem_Translation_Method()( *(ch.chrom), pr );
                        ch.fit = Evaluation_Method()(pr);
                        ch.modification_after_last_evaluation_flag = false;
                        current_population.push_back( ch );
                    }
                }
                else
                {
                    //shuffle the resize
                    std::random_shuffle( current_population.begin(), current_population.end() );
                    current_population.resize( population_per_generation );
                }
    
                std::sort( current_population.begin(), current_population.end() );

                //total_fit = 0;
                //for( auto& ch : current_population )
                //    total_fit += ch.fit;
                //std::cout.precision(20);
                //std::cout << total_fit << "\n";

                //if timeout, return output elite
                auto& t_manager = feng::singleton<time_manager>::instance();
                if ( t_manager.is_timeout() )
                {
                    std::cerr << "\nthe residual for the best individual is " << best_fitness;

                    elite_out.close();
                    elite_fit.close();

                    return best_one;
                    //return pr;
                }

                //std::cerr << "\nElite of " << debug_counter++ << " is \n" << pr;
                //std::cerr << "\nElite of " << debug_counter << " is \n" << *(elite_one.chrom);


            //select 
                auto& xpm = feng::singleton<xover_probability_manager<>>::instance();
                std::size_t const selection_number = xpm(population_per_generation);
                //selected_population_for_xover_father.resize(selection_number);
                //selected_population_for_xover_mother.resize(selection_number);
                selected_population_for_xover_father.clear();
                selected_population_for_xover_mother.clear();

                //std::cerr << "\nselecting............." << selection_number;

                auto& xs_manager = feng::singleton<xover_selection_manager>::instance();
                for ( std::size_t i = 0; i != selection_number; ++i )
                {
                    //selected_population_for_xover_father[i] = current_population[xs_manager()];
                    //selected_population_for_xover_mother[i] = current_population[xs_manager()];
                    
                    const std::size_t select1 = xs_manager();
                    const std::size_t select2 = xs_manager();

                    //selected_population_for_xover_father.push_back(current_population[xs_manager()]);
                    //selected_population_for_xover_mother.push_back(current_population[xs_manager()]);
                    selected_population_for_xover_father.push_back(current_population[select1]);
                    selected_population_for_xover_mother.push_back(current_population[select2]);
                }

                //std::cerr << "\nselectedted";
            //xover
#if 0
                //not working as selected population for xover are pure reference of current populaiton
                for ( std::size_t i = 0; i != selection_number; ++i )
                    feng::for_each( selected_population_for_xover_father[i].begin(), selected_population_for_xover_father[i].end(), 
                                    selected_population_for_xover_mother[i].begin(), Chromosome_Xover_Method() );

                current_population.reserve( population_per_generation + selection_number + selection_number );
                //append newly generated genes into current population
                current_population.insert( current_population.begin()+population_per_generation, selected_population_for_xover_father.begin(), selected_population_for_xover_father.end() );
                current_population.insert( current_population.begin()+population_per_generation+selection_number, selected_population_for_xover_mother.begin(), selected_population_for_xover_mother.end() );
#endif
                //std::cerr << "\nxover.............";
                
                for ( std::size_t i = 0; i != selection_number; ++i )
                {
                    chromosome_type son( ch_args... ); 
                    chromosome_type daughter( ch_args... ); 

                    feng::for_each( selected_population_for_xover_father[i].begin(), selected_population_for_xover_father[i].end(), selected_population_for_xover_mother[i].begin(), son.begin(), daughter.begin(), Chromosome_Xover_Method() );

                    current_population.push_back( son ); 
                    current_population.push_back( daughter ); 
                }
                //mark new generation not evaluated
                feng::for_each( current_population.begin()+population_per_generation, current_population.end(), [](chromosome_type& ch) { ch.modification_after_last_evaluation_flag = true; } );

                //std::cerr << "\nxovered";

            //mutate
                //std::cerr << "\nmutating";

                auto& mpm = feng::singleton<mutation_probability_manager<>>::instance();
                mpm.reset();
                for ( auto& chromosome : current_population )
                    if ( mpm.should_mutate_current_gene() )
                    {
                        for ( auto& gene : chromosome )
                            Chromosome_Mutation_Method()(gene);
                        //makr muated ones as not evaluated
                        chromosome.modification_after_last_evaluation_flag = true;
                    }

                //std::cerr << "\nmutated";

            //goto evaluate



            }

            assert( !"Should never reach here!!" );

            return Problem_Representation_Type();
                


        }

    };//struct genetic_algorithm 

}//namespace ga

#endif//_GENETIC_ALGORITHM_HPP_INCLUDED_SFOJ3W94889YKLSFDJHXVCMSAJKL4389SFDJKSFAJKLSFAFIOOSFIDWREOFF

