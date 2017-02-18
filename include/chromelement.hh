//
// Abstract building block of Chromosome.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_CHROMELEMENT_H_
#define _DORSALFIN_CHROMELEMENT_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class ChromosomeElement
    /// \brief Abstract element of a chromosome.
    ///
    /// In \c fluke chromosomes consist of \em atoms, called chromosome 
    /// elements. 
    class ChromosomeElement : public CachedElement, public XmlWriter {
        public:
        /// Virtual destructor.
        virtual ~ChromosomeElement() {};
        /// Signature of clone method.
        virtual ChromosomeElement* clone() const = 0;
        /// Signature of copy method.
        virtual void copy( const ChromosomeElement & ) = 0;
        
        /// Deactivate the element.
        void inactivate();
        /// Activate the element.
        void activate();
        /// Is the element \em activated?
        bool isActive() const;
        
        /// Write to an output stream.
        void write( std::ostream & ) const;
        /// Write to an xml stream
        virtual void writeXml( std::ostream & ) const = 0;
        /// write a string representation
        //virtual std::string asString() const = 0;

        protected:
        /// Hidden constructor.
        ChromosomeElement() : CachedElement(), active_( true ) {};
        /// Hidden copy constructor.
        explicit ChromosomeElement( const ChromosomeElement &ce ) 
            : CachedElement(), active_( ce.active_ ) {};

        private:
        /// Flag signaling if the element is active.
        bool active_;
    };

    inline void ChromosomeElement::inactivate()
    { active_ = false; }

    inline void ChromosomeElement::activate()
    { active_ = true; }

    inline bool ChromosomeElement::isActive() const
    { return active_; }

    inline void ChromosomeElement::write( std::ostream &os ) const 
    { writeXml( os ); }
}
#endif

