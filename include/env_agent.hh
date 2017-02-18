//
// Interface of an simple agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_ENVAGENT_H_
#define _DORSALFIN_ENVAGENT_H_

#include "defs.hh"
#include "genreg_agent.hh"
#include "genome.hh"
#include "network.hh"

namespace dorsalfin {

    /// \class EnvAgent
    /// \brief Simple agent with a network.
    class EnvAgent : public GenRegAgent {
        public:
        typedef std::vector< boost::dynamic_bitset<> >::iterator vb_iter;

        public:
        /// Default ctor.
        EnvAgent();
        /// Ctor.
        EnvAgent( uint tag, Genome *g );
        /// Copy ctor.
        EnvAgent( const EnvAgent & );
        /// Destructor.
        virtual ~EnvAgent();
        /// Cloning an agent.
        virtual Agent* clone() const;
        /// Copy the agent into \c this
        virtual void copy( const Agent & );

        /// Perform one timestep. It consists of checking whether the
        /// agent dies or survives another timestep.
        virtual void step( Population & );
        /// Reproduce. The sibling has its age reset to zero.
        virtual Agent* sibling();
        /// Evaluate network and calculate new attractor
        virtual void evaluate( const Environment &, const location & );

        /// Overriding score method
        virtual double score() const;
        
        /// Reproduction?
        bool reproduced() const;
        /// What was the input sequence?
        boost::dynamic_bitset<> input() const;
        /// Chopping of any leading zeroes/ones
        boost::dynamic_bitset<> output() const;
        /// How many bits were actually correct (w/ max at nr of input genes)
        int nrPredictedBits() const;
        
        /// Write a text representation of the agent to an output stream.
        virtual void writeXml( std::ostream & ) const;
    
        private:
        // produced offspring
        bool reproduced_;
    };
    
    inline bool EnvAgent::reproduced() const
    { return reproduced_; }
    
    inline int EnvAgent::nrPredictedBits() const 
    { return distance_; }

    inline boost::dynamic_bitset<> dorsalfin::EnvAgent::input() const 
    { return network_->input(); }

    inline boost::dynamic_bitset<> dorsalfin::EnvAgent::output() const 
    { return network_->output(); }
}
#endif

