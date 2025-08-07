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
	bool operator() (const Box* a, const Box* b) const
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

	BoxQueue(void)
	{
	}

	virtual void push(Box* b) = 0;

	/** Pops off the item at the front of the queue and returns it */
    virtual Box* extract() = 0;

	virtual bool empty() = 0;

	virtual int size() = 0;

	~BoxQueue(void)
	{
	}
};

/*! \class seqQueue
 *  \brief "Normal" (not randomized) priority queue of `Box` nodes
 */ 
class seqQueue : public BoxQueue
{
private:
	priority_queue<Box*, vector<Box*>, PQCmp> PQ;
public:
	void push(Box* b)
	{
		PQ.push(b);
	}

    /** Pops off the item at the front of the queue and returns it */
	Box* extract()
	{
		Box* r = PQ.top();
		PQ.pop();
		return r;
	}

	bool empty()
	{
		return PQ.empty();
	}

	int size()
	{
		return PQ.size();
	}
};

/*! \class randQueue
 *  \brief List of `Box` nodes with randomization
 *
 *  A list of `Box` nodes, where list extraction
 *  returns a random node.
 */ 
class randQueue : public BoxQueue
{
private:
	list<Box*> L;
	int Qseed;

public:
	randQueue(int s): Qseed(s) {
		//srand( time(0) );
		srand( Qseed ); 
	}

	void push(Box* b)
	{
		L.push_back(b);
	}

	/** Removes and returns a random item from the list */
    Box* extract()
	{
		int i = rand() % L.size();
		list<Box*>::iterator iter = L.begin();
		advance(iter, i);
		Box* r = *iter;
		L.erase(iter);
		return r;
	}

	bool empty()
	{
		return L.empty();
	}

	int size()
	{
		return L.size();
	}

};


} //end namespace cusg
