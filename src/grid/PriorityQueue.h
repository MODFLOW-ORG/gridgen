/*! \file PriorityQueue.h
 *
 *  \brief Priority queue of `Box` nodes
 *  \author Jyh-Ming Lien
 *  \author <a href="http://masc.cs.gmu.edu">MASC Group</a>, 
 *          George Mason University
 *  \bug No known bugs.
 */

//$Id: $
//
#pragma once

#include "Box.h"
#include <queue>
#include <vector>
#include <list>
#include <time.h>
#include <stdlib.h>
#include <iterator>
#include <math.h>

using namespace std;

namespace cusg
{

/*! \class PQCmp 
 *  \brief `Box` comparer
 *
 * Used to compare `Box` nodes in a priority queue.
 * Boxes are compared by depth, then by time of creation.
 */
class PQCmp
{
public:
    bool operator()(const Box * a, const Box * b) const
    {
        //use depth for now
        if (a->depth > b->depth)
        {
            return true;
        }
        //if same depth, expand box created earlier first
        else if (a->depth == b->depth)
        {
            return a->priority > b->priority;
        }
        return false;
    }
};

/*! \class BoxQueue 
 *  \brief Dummy base class for `Box` node collection classes */
class BoxQueue
{
private:
public:
    BoxQueue(void) {}

    virtual void push(Box * b) = 0;

    /** Pops off the item at the front of the queue and returns it */
    virtual Box * extract() = 0;

    virtual bool empty() = 0;

    virtual int size() = 0;

    ~BoxQueue(void) {}
};

/*! \class SeqQueue
 *  \brief "Normal" (not randomized) priority queue of `Box` nodes
 */
class SeqQueue : public BoxQueue
{
private:
    priority_queue<Box *, vector<Box *>, PQCmp> m_pq;

public:
    void push(Box * b) override { m_pq.push(b); }

    /** Pops off the item at the front of the queue and returns it */
    Box * extract() override
    {
        Box * r = m_pq.top();
        m_pq.pop();
        return r;
    }

    bool empty() override { return m_pq.empty(); }

    int size() override { return m_pq.size(); }
};

/*! \class RandQueue
 *  \brief List of `Box` nodes with randomization
 *
 *  A list of `Box` nodes, where list extraction
 *  returns a random node.
 */
class RandQueue : public BoxQueue
{
private:
    list<Box *> m_l;
    int m_qseed;

public:
    RandQueue(int s) : m_qseed(s) { srand(m_qseed); }

    void push(Box * b) override { m_l.push_back(b); }

    /** Removes and returns a random item from the list */
    Box * extract() override
    {
        int i = rand() % m_l.size();
        list<Box *>::iterator iter = m_l.begin();
        advance(iter, i);
        Box * r = *iter;
        m_l.erase(iter);
        return r;
    }

    bool empty() override { return m_l.empty(); }

    int size() override { return m_l.size(); }
};

} //end namespace cusg
