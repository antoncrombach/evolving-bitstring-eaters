//
// Caching objects.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_POOL_H_
#define _DORSALFIN_POOL_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class CachedElement
    /// \brief Interface for cache-able object
    class CachedElement {
        public:
        /// Virtual dtor.
        virtual ~CachedElement() {};
        
        /// Signature of return to pool.
        virtual bool toPool() = 0;
        
        protected:
        /// Hidden ctor.
        CachedElement() {};
        /// Hidden copy ctor.
        explicit CachedElement( const CachedElement &p ) {};
    };


    /// \class ObjectCache
    /// \brief The class caches small objects
    ///
    /// The creation and deletion of thousands (and more) small objects
    /// is one of the bottlenecks of my current implementation. My hope
    /// is that pooling or caching such small instances speeds up te simulation
    /// considerably.
    ///
    /// Note: uses pointers instead of entire instances.
    /// Question: how to use a factory with this pool?
    template< class T >
    class ObjectCache {
        public:
        /// Dtor
        ~ObjectCache();
        
        /// Get an object. If pool is empty, a new object is created and 
        /// returned.
        T* borrowObject();
        /// Return an idle instance to the pool. If the maximum of idle
        /// instances is reached, just return false.
        bool returnObject( T* );
        
        /// Add an idle instance to the pool.
        void addObject();
        /// Clear the pool. All instances are destroyed.
        void clear();
        
        /// Return number of active instances. These are not under control of
        /// the pool, so provide an upperbound.
        uint getNrActive() const;
        /// Return number of idle instances.
        uint getNrIdle() const;
        
        public:
        /// Get the pool.
        static ObjectCache< T >* instance();

        protected:
        /// Constructor.
        ObjectCache();
        /// Constructor with maximum number of idle instances.
        ObjectCache( uint );
        
        protected:
        /// Number of active instances (outside the pool).
        uint nr_active_;
        /// Cap number of idle instances. 
        uint max_idle_;
        /// The pool itself.
        std::stack< T* > pool_;
        
        protected:
        static ObjectCache< T > *instance_;
    };
    
    // Note: using magic number
    template< class T >
    ObjectCache< T >::ObjectCache()
        : nr_active_( 0 ),  max_idle_( 16 ), pool_() {}
    
    template< class T >
    ObjectCache< T >::ObjectCache( uint maxIdle ) 
        : nr_active_( 0 ),  max_idle_( maxIdle ), pool_() {}

    template< class T >
    ObjectCache< T >::~ObjectCache() {
        clear();
    }
        
    template< class T > T*
    ObjectCache< T >::borrowObject() {
        ++nr_active_;
        if( pool_.empty() ) {
            addObject();
        }
        T* result = pool_.top();
        pool_.pop();
        return result;
    }
    
    template< class T > bool
    ObjectCache< T >::returnObject( T* obj ) {
        --nr_active_;
        if( pool_.size() == max_idle_ ) {
            //delete obj;
            return false;
        } else {
            pool_.push( obj );
            return true;
        }
    }
    
    template< class T > void
    ObjectCache< T >::addObject() {
        T *obj = new T();
        pool_.push( obj );
    }
    
    template< class T > void 
    ObjectCache< T >::clear() {
        while( !pool_.empty() ) {
            delete pool_.top();
            pool_.pop();
        }
    }

    template< class T > uint ObjectCache< T >::getNrActive() const
    { return nr_active_; }
    
    template< class T > uint ObjectCache< T >::getNrIdle() const
    { return max_idle_; }
    
    // Note: using magic number
    template< class T > ObjectCache< T >* 
    ObjectCache< T >::instance() {
        if( instance_ == 0 ) {
            instance_ = new ObjectCache< T >( 5 * 16384 );
        }
        return instance_;
    }
    
    /// General definition
    template< class T > ObjectCache< T >* ObjectCache< T >::instance_ = 0;
    
    /// Extra function related to the use of \c ObjectCache
    ///
    /// Return to pool (use ONLY for classes derived from CachedElement).
    template< class Container, class For > For
    smart_return( Container &x, For first ) {
        if( !( **first ).toPool() ) delete *first;
        return x.erase( first );
    }

    /// Extra function related to the use of \c ObjectCache
    ///
    /// Range return to pool (use ONLY for classes derived from CachedElement).
    /// See also smart_return( Container, For )
    template< class Container, class For > For
    smart_return( Container &x, For first, For last ) {
        For i = first;
        while( i != last ) {
            if( !( **i ).toPool() ) delete *i;
            ++i;
        }
        return x.erase( first, last );
    }

}
#endif
