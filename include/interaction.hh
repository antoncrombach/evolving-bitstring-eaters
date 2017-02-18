//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_INTERACTION_H_
#define _DORSALFIN_INTERACTION_H_

#include "defs.hh"
#include "pool.hh"
#include "chromelement.hh"

namespace dorsalfin {

    /// \class Interaction
    /// \brief Interactions between genes.
    ///
    /// Binding sites in the upstream region of a genome act as the 'glue' in
    /// interactions between genes. Transcription Factors bind to these sites
    /// activating or inhibiting the downstream gene.
    class Interaction : public ChromosomeElement {
        public:
        /// Constructor.
        Interaction();
        /// Constructor with arguments
        Interaction( int, int );
        /// Copy constructor.
        explicit Interaction( const Interaction & );
        /// Destructor.
        virtual ~Interaction() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another Interaction into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return Interaction to pool
        virtual bool toPool();
        
        /// Interaction is activating (+1) or inhiting (-1)
        void weight( int );
        /// Interaction weight (+-1)
        int weight() const;
        /// Toggle the interaction value
        void toggle();
        
        /// The tag signals which genes bind here 
        /// (there may be multiple copies of the genes)
        void tag( int );
        /// The tag 
        int tag() const;
        
        /// Write a html representation.
        virtual void writeXml( std::ostream & ) const;
        /// Write a string representation
        virtual std::string asString() const;
        
        private:
        int weight_;
        int tag_;
    };

    inline void Interaction::weight( int w )
    { weight_ = w; }
    
    inline int Interaction::weight() const
    { return weight_; }
    
    inline void Interaction::toggle()
    { weight_ = -weight_; }
    
    inline void Interaction::tag( int t )
    { tag_ = t; }
    
    inline int Interaction::tag() const
    { return tag_; }
    
    inline std::string Interaction::asString() const
    { return "ia"; }
}
#endif

