//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVertex - A node in a scaffold graph. 
//
#include "ScaffoldVertex.h"

//
ScaffoldVertex::ScaffoldVertex(VertexID id, size_t seqLen) : m_id(id), m_seqLen(seqLen), m_classification(SVC_UNKNOWN), m_color(GC_WHITE)
{

}

//
ScaffoldVertex::~ScaffoldVertex()
{
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
          iter != m_edges.end(); ++iter)
    {
        delete *iter;
        *iter = 0;
    }
}

//
void ScaffoldVertex::addEdge(ScaffoldEdge* pEdge)
{
    m_edges.push_back(pEdge);
}

//
ScaffoldEdge* ScaffoldVertex::findEdgeTo(VertexID id, ScaffoldEdgeType type) const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id && (*iter)->getType() == type)
            return *iter;
    }

    return NULL;
}

//
ScaffoldEdge* ScaffoldVertex::findEdgeTo(VertexID id, EdgeDir dir, EdgeComp comp) const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id && (*iter)->getDir() == dir && (*iter)->getComp() == comp)
            return *iter;
    }

    return NULL;
}

//
ScaffoldEdgePtrVector ScaffoldVertex::getEdges()
{
    return m_edges;
}

//
ScaffoldEdgePtrVector ScaffoldVertex::getEdges(EdgeDir dir)
{
    ScaffoldEdgePtrVector out;
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDir() == dir)
            out.push_back(*iter);
    }
    return out;
}

void ScaffoldVertex::deleteEdges()
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        delete *iter;
        *iter = NULL;
        ++iter;
    }
    m_edges.clear();
}

void ScaffoldVertex::deleteEdgesAndTwins()
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        ScaffoldEdge* pEdge = *iter;
        pEdge->getEnd()->deleteEdge(pEdge->getTwin());
        delete pEdge;
        *iter = NULL;
        ++iter;
    }
    m_edges.clear();
}

void ScaffoldVertex::deleteEdgesAndTwins(EdgeDir dir)
{
    ScaffoldEdgePtrVector out;

    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        ScaffoldEdge* pEdge = *iter;

        // We have to clean up the twin of each deleted edge. If the twin
        // points to this vertex (self-edge) we can't clean it up with a call to deleteEdge
        // or else our iterator will be invalidated. We handle it here by deleting
        // each edge if it matches the direction or it is the twin of an edge that
        // matches the direction. In the inner body we skip the standard twin cleanup
        // for self-edges
        if(pEdge->getDir() == dir || (pEdge->getEnd() == this && pEdge->getTwin()->getDir() == dir))
        {
            // If the twin resides on this node and points in the same direction
            // it will be cleaned up in a later cycle so do not delete it
            if(pEdge->getEnd() != this)
            {
                pEdge->getEnd()->deleteEdge(pEdge->getTwin());
            }
            delete pEdge;
            *iter = NULL;
        }
        else
        {
            out.push_back(*iter);
        }
        ++iter;
    }
    
    m_edges = out;
}

// Remove an edge from the edge vector and delete it
void ScaffoldVertex::deleteEdge(ScaffoldEdge* pEdge)
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if(*iter == pEdge)
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);
    delete pEdge;
}

// Remove an edge from the edge vector and delete it and its twin
void ScaffoldVertex::deleteEdgeAndTwin(ScaffoldEdge* pEdge)
{
    assert(pEdge != pEdge->getTwin());
    pEdge->getEnd()->deleteEdge(pEdge->getTwin());
    deleteEdge(pEdge);
}

//
void ScaffoldVertex::setAStatistic(double v)
{
    m_AStatistic = v;
}

//
void ScaffoldVertex::setClassification(ScaffoldVertexClassification classification)
{
    m_classification = classification;
}

//
void ScaffoldVertex::setColor(GraphColor c)
{
    m_color = c;
}

//
VertexID ScaffoldVertex::getID() const
{
    return m_id;
}

//
bool ScaffoldVertex::isRepeat() const
{
    return m_classification == SVC_REPEAT;
}
//
size_t ScaffoldVertex::getNumEdges() const
{
    return m_edges.size();
}

//
size_t ScaffoldVertex::getSeqLen() const
{
    return m_seqLen;
}

//
double ScaffoldVertex::getAStatistic() const
{
    return m_AStatistic;
}

//
ScaffoldVertexClassification ScaffoldVertex::getClassification() const
{
    return m_classification;
}

//
GraphColor ScaffoldVertex::getColor() const
{
    return m_color;
}

//
std::string ScaffoldVertex::getColorString() const
{
    switch(m_classification)
    {
        case SVC_UNKNOWN:
            return "gray";
        case SVC_UNIQUE:
            return "white";
        case SVC_REPEAT:
            return "red";
        default:
            return "white";
    }
}

//
void ScaffoldVertex::writeDot(std::ostream* pWriter) const
{
   VertexID id = getID();
   *pWriter << "\"" << id << "\" [ label =\"" << id << "," << getSeqLen() << "\" ";
   *pWriter << "style=\"filled\" fillcolor=\"" << getColorString() << "\" ";
   *pWriter << "];\n";
   writeEdgesDot(pWriter);
}

//
void ScaffoldVertex::writeEdgesDot(std::ostream* pWriter) const
{
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        *pWriter << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
        std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
        *pWriter << "\" [";
        *pWriter << "label=\"" << (*iter)->getDistance() << "\" ";
        *pWriter << "color=\"" << color << "\" ";
        *pWriter << "];";
        *pWriter << "\n";
    }
}