/*! \file CubeQueue.h
 *
 *  \brief Priority queue of `Cube` nodes
 *  \author Jyh-Ming Lien
 *  \author <a href="http://masc.cs.gmu.edu">MASC Group</a>, 
 *          George Mason University
 *  \bug No known bugs.
 */

//$Id: $
//
#pragma once

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <list>
#include <queue>
#include <vector>
#include "Cube.h"

using namespace std;

namespace cusg
{

/*! \class CubeQueueComparer
 *  \brief `Cube` comparer
 *
 * Used to compare `Cube` nodes in a priority queue.
 * Cubes are compared by depth, then by time of creation.
 */ 
class CubeQueueComparer
{
public:
	bool operator() (const Cube* a, const Cube* b) const
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

/*! \class CubeQueue 
 *  \brief Dummy base class for `Cube` node collection classes */
class CubeQueue
{
private:

public:

	CubeQueue(void)
	{
	}

	virtual void push(Cube* b) = 0;

	/** Pops off the item at the front of the queue and returns it */
    virtual Cube* extract() = 0;

	virtual bool empty() = 0;

	virtual int size() = 0;

	~CubeQueue(void)
	{
	}
};

/*! \class SequentialCubeQueue
 *  \brief "Normal" (not randomized) priority queue of `Cube` nodes
 */ 
class SequentialCubeQueue : public CubeQueue
{
private:
	priority_queue<Cube*, vector<Cube*>, CubeQueueComparer> PQ;
public:
	void push(Cube* b)
	{
		PQ.push(b);
	}

    /** Pops off the item at the front of the queue and returns it */
	Cube* extract()
	{
		Cube* r = PQ.top();
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

/*! \class RandomizedCubeQueue
 *  \brief List of `Cube` nodes with randomization
 *
 *  A list of `Cube` nodes, where list extraction
 *  returns a random node.
 */ 
class RandomizedCubeQueue : public CubeQueue
{
private:
	list<Cube*> L;
	int Qseed;

public:
	RandomizedCubeQueue(int s): Qseed(s) {
		srand(Qseed); 
	}

	void push(Cube* b)
	{
		L.push_back(b);
	}

	/** Removes and returns a random item from the list */
    Cube* extract()
	{
		int i = rand() % L.size();
		list<Cube*>::iterator iter = L.begin();
		advance(iter, i);
		Cube* r = *iter;
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
