#pragma once
#ifndef CUSG_QUADTREE_SHP_H
#define CUSG_QUADTREE_SHP_H

#include "QuadTree3D.h"
#include "shapelib/shapefil.h"
#include <string>
#include <algorithm>
#include <vector>
using namespace std;
using namespace cusg;

//write a quadtree to a shape file
void writeQuadtree2Shp(QuadTree3D* qtree, string shpName, string feature_type, bool bConcide=true);

//write a modflow grid to a shape file
void writeModflowGrid2Shp(ModflowGrid* modflow, string gridShpName, string feature_type,bool without_inactive = true, bool one_based_numbering = true, bool bConcide=true);



////////////////////////////////////////////////////////////////////////////////
// Implementation of sorting return index
////////////////////////////////////////////////////////////////////////////////
// Input:
//   unsorted  unsorted vector
// Output:
//   index_map  an index map such that sorted[i] = unsorted[index_map[i]]
template <class T>
void sort(
    std::vector<T> &unsorted,
    std::vector<size_t> &index_map);

// Comparison struct used by sort
template<class T> struct index_cmp 
{
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  { 
    return arr[a] < arr[b];
  }
  const T arr;
};

template <class T>
void sort(
  std::vector<T> & unsorted,
  std::vector<size_t> & index_map)
{
  // Original unsorted index map
  index_map.resize(unsorted.size());
  for(size_t i=0;i<unsorted.size();i++)
  {
    index_map[i] = i;
  }
  // Sort the index map, using unsorted for comparison
  sort(
    index_map.begin(), 
    index_map.end(), 
    index_cmp< std::vector<T> & >(unsorted));
}


#endif
