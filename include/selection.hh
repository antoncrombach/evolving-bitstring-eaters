//
// Interface of the population selection schemes.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_SELECTIONSCHEME_H_
#define _DORSALFIN_SELECTIONSCHEME_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class SelectionScheme
    /// \brief An agent is selected for reproduction on basis of its score
    ///
    /// Given a list of agents and a list of scores (and assuming both are 
    /// ordered equivallently), one agent is selected.
    class SelectionScheme {
        public:
            /// Destructor
            virtual ~SelectionScheme() {};
            /// Cloning is very handy
            virtual SelectionScheme* clone() const = 0;
            
            /// Signature of the select method
            virtual Agent* select( 
                    std::vector< Agent* > &, std::vector< double > & ) = 0;

        protected:
            /// Constructor
            SelectionScheme() {};
    };

    /// \class RandSelection
    /// \brief Randomly pick an agent
    ///
    /// The easiest way of selecting an agent is ignoring the scores and just
    /// randomly pick one of them.
    class RandSelection : public SelectionScheme {
        public:
            /// Constructor
            RandSelection() : SelectionScheme() {};
            /// Destructor
            virtual ~RandSelection() {};

            virtual SelectionScheme* clone() const;

            /// Select an agent (randomly)
            virtual Agent* select( 
                    std::vector< Agent* > &, std::vector< double > & );
    };

    inline SelectionScheme* RandSelection::clone() const
    { return new RandSelection(); }
    
    /// \class ProbalisticSelection
    /// \brief Pick an agent linearly based on its score
    ///
    /// An agent is selected on basis of its score. Agents with a higher score
    /// have a higher probability of being picked, yet all agents have a chance
    /// of being selected. 
    class ProbalisticSelection : public SelectionScheme {
        public:
            /// Constructor
            ProbalisticSelection() : SelectionScheme() {};
            /// Destructor
            virtual ~ProbalisticSelection() {};

            virtual SelectionScheme* clone() const;

            /// Select an agent on its score
            virtual Agent* select( 
                    std::vector< Agent* > &, std::vector< double > & );
    };
            
    inline SelectionScheme* ProbalisticSelection::clone() const
    { return new ProbalisticSelection(); }
}
#endif

