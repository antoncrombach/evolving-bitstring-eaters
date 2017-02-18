//
// Abstract class of Gene regions
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_GENE_H_
#define _DORSALFIN_GENE_H_

#include "defs.hh"
#include "pool.hh"
#include "chromelement.hh"
#include "ref_node.hh"

namespace dorsalfin {

    /// \class Gene
    /// \brief Representation of a gene
    class Gene : public ChromosomeElement {
        public:
        /// Constructor.
        Gene();
        /// Constructor with tag
        Gene( int );
        /// Constructor with tag, state and threshold
        Gene( int, int, int );
        /// Copy constructor.
        explicit Gene( const Gene & );
        /// Virtual destructor.
        virtual ~Gene() {}
        // Cloning.
        virtual ChromosomeElement* clone() const;
        /// Copy method.
        virtual void copy( const ChromosomeElement & );
        /// Return to pool method.
        virtual bool toPool();

        /// Mutate threshold
        void mutateThreshold();
        /// Mutate state
        void mutateState();
        /// Import data from reference network
        void import( const RefNode & );
        
        /// Set tag
        void tag( int );
        /// Get tag
        int tag() const;
        /// Set state
        void state( int );
        /// Get state
        int state() const;
        /// Set threshold
        void threshold( int );
        /// Get threshold of expression
        int threshold() const;

        /// Write a html representation.
        virtual void writeXml( std::ostream & ) const;
        /// Write a string representation
        virtual std::string asString() const;
        /// Write a label for graphs.
        virtual std::string asLabel() const;
        
        protected:
        /// Identification and gives target bsite type
        int tag_;
        /// On or off
        int state_;
        /// Regulation threshold, has a value {-2,-1,0,1,2}]
        /// Note: these values have been hardcoded.......
        int threshold_;
    };

    inline void Gene::tag( int t )
    { tag_ = t; }

    inline int Gene::tag() const
    { return tag_; }
    
    inline void Gene::state( int t )
    { state_ = t; }
    
    inline int Gene::state() const
    { return state_; }
    
    inline void Gene::threshold( int t )
    { threshold_ = t; }

    inline int Gene::threshold() const 
    { return threshold_; }
    
    inline std::string Gene::asLabel() const
    { return asString(); }
    
    inline std::string Gene::asString() const
    { return boost::lexical_cast< std::string >( tag_ ); }
    
    
}
#endif
