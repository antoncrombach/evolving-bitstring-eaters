//
// Reference network needs special node
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_REFERENCE_NODE_H_
#define _DORSALFIN_REFERENCE_NODE_H_

namespace dorsalfin {
    
    /// \class RefNode
    /// \brief Auxilary class for reading in reference network
    class RefNode {
        public:
        RefNode() : node_id( 0 ), threshold( 0 ), 
            transfac( false ), input( false ) {}
        RefNode( const RefNode &rf ) : node_id( rf.node_id ), 
            threshold( rf.threshold ), transfac( rf.transfac ), 
            input( rf.input ) {}
        ~RefNode() {}
        
        public:
        int node_id;
        int threshold;
        bool transfac;
        bool input;
    };
}
#endif
