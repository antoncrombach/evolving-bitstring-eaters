//
// Interface of the population scaling schemes.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_SCALINGSCHEME_H_
#define _DORSALFIN_SCALINGSCHEME_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class ScalingScheme
    /// \brief The agents score is scaled to the interval \f$[0,1]\f$
    ///
    /// Reproduction is based on the score of agents. As in genetic algorithms
    /// first the fitness score is scaled.
    class ScalingScheme {
        public:
            /// Destructor
            virtual ~ScalingScheme() {};
            /// Cloning is very handy
            virtual ScalingScheme* clone() const = 0;
    
            /// Signature of the scaling method
            virtual void scale( std::vector< double > & ) = 0;

        protected:
            /// Constructor
            ScalingScheme() {};
    };

    /// \class NoScaling
    /// \brief Identity scaling
    ///
    /// The simplest form of scaling is no scaling at all.
    class NoScaling : public ScalingScheme {
        public:
            /// Constructor
            NoScaling() : ScalingScheme() {};
            /// Destructor
            virtual ~NoScaling() {};

            virtual ScalingScheme* clone() const;

            /// Scale the scores, but not really ;)
            virtual void scale( std::vector< double > & );
    };

    inline ScalingScheme* NoScaling::clone() const
    { return new NoScaling(); }
    
    /// \class LinearScaling
    /// \brief Linearly scale the scores from \f$[0, \infty)\f$ to \f$[0,1]\f$
    ///
    /// The scores are scaled according to the formula 
    /// \f$\frac{score}{score_{max}}\f$. In case all scores are equal, a base
    /// score is given to each score. This is the scaled score.
    class LinearScaling : public ScalingScheme {
        public:
            /// Constructor
            LinearScaling() : ScalingScheme(), base_score_( 0.0 ) {};
            /// Constructor. The double \c f is usually the birth rate of agents
            LinearScaling( double f ) : ScalingScheme(), base_score_( f ) {};
            /// Destructor
            virtual ~LinearScaling() {};

            virtual ScalingScheme* clone() const;
            
            /// Scale the scores linearly
            virtual void scale( std::vector< double > & );
        private:
            double base_score_;
    };

    inline ScalingScheme* LinearScaling::clone() const
    { return new LinearScaling( base_score_ ); }

    /// \class PowerScaling
    /// \brief All scores are raised to a specified power (default 2).
    class PowerScaling : public ScalingScheme {
        public:
            /// Constructor
            PowerScaling() 
                : ScalingScheme(), base_score_( 0.0 ), power_( 2.0 ) {};
            /// Constructor. The double \c f is usually the birth rate of agents
            PowerScaling( double f, double p ) 
                : ScalingScheme(), base_score_( f ), power_( p ) {};
            /// Destructor
            virtual ~PowerScaling() {};

            virtual ScalingScheme* clone() const;
            
            /// Scale the scores linearly
            virtual void scale( std::vector< double > & );
        private:
            double base_score_;
            double power_;
    };

    inline ScalingScheme* PowerScaling::clone() const
    { return new PowerScaling( base_score_, power_ ); }

    /// \class ExponentialScaling
    /// \brief All scores are exponentially scaled (with a coefficient).
    class ExponentialScaling : public ScalingScheme {
        public:
            /// Constructor
            ExponentialScaling() 
                : ScalingScheme(), base_score_( 0.0 ), power_( 1.0 ) {};
            /// Constructor. The double \c f is usually the birth rate of agents
            ExponentialScaling( double f, double p ) 
                : ScalingScheme(), base_score_( f ), power_( p ) {};
            /// Destructor
            virtual ~ExponentialScaling() {};

            virtual ScalingScheme* clone() const;
            
            /// Scale the scores linearly
            virtual void scale( std::vector< double > & );
        private:
            double base_score_;
            double power_;
    };

    inline ScalingScheme* ExponentialScaling::clone() const
    { return new ExponentialScaling( base_score_, power_ ); }
}
#endif

