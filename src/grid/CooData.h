
/*! \file CooData.h

\brief A structure to hold a sparse matrix in Coo format

\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/


#pragma once

#include <cassert>
#include "CsrData.h"
using namespace std;

namespace cusg
{

/*! \brief Sparse Matrix in Coo format

    This object contains all of the unstructured grid information
    necessary for an unstructured grid model.  The data are stored
    in the Coo format.
 */
struct CooData
{
	///
	CooData(){ } 

    ///Constructor
    CooData(int nodes, int ncon, int * coorow, int * coocol, int * indptr, double * area, float * cl1, float * cl2, float * fahl, int * fldr)
    {
        this->nodes = nodes;
        this->ncon = ncon;
        this->areas = area;
        this->coorow = coorow;
        this->coocol = coocol;
        this->indptr = indptr;
        this->cl1 = cl1;
        this->cl2 = cl2;
        this->fahl = fahl;
        this->fldr = fldr;
    }

    /// Creates a matrix in Coo format, given a sparse matrix in CSR format.
    CooData(CSRData csr)
    {
        int * coorow = new int[csr.nja];
        int * coocol = new int[csr.nja];
        int *  indptr = new int[csr.nodes+1];

        assert(coorow && coocol && indptr);

        int ipos = 0;
        for( int n=0; n<csr.nodes; n++)
        {
            indptr[n] = ipos;

            for(int j=csr.ia[n];j<csr.ia[n + 1];j++)
            {
                int m = csr.ja[j-1]; //TODO:  -1 b/c csr is 1-based
                coorow[ipos] = n;
                coocol[ipos] = m-1; //TODO:  -1 b/c csr is 1-based
				//cout<<ipos<<" ";
                ipos += 1;
            }
        }

        this->nodes = csr.nodes;
        this->ncon = csr.nja;
        this->coorow = coorow;
        this->coocol = coocol;
        this->indptr = indptr;
        this->areas = csr.area;
        this->cl1 = csr.cl1;
        this->cl2 = csr.cl2;
        this->fahl = csr.fahl;
        this->fldr = csr.fldr;
    }

    ///Remove diagonal entries from the arrays and reduce size.
    void remove_diagonal()
    {
        //count nondiagonal entries and reallocate
        int ncon = 0;
        for( int i=0;i<this->ncon;i++)
        {
            if( this->coorow[i] != this->coocol[i]) ncon += 1;
        }

        int * coorow = new int[ncon];
        int * coocol = new int[ncon];
        float * cl1 = new float[ncon];
        float * cl2 = new float[ncon];
        float * fahl = new float[ncon];
        int * fldr = new int[ncon];

        assert(coorow && coocol && cl1 && cl2 && fahl && fldr);

        //reassign
        int inewpos = 0;
        int lastn = -1;
        for( int i=0;i<this->ncon;i++)
        {
            if( this->coorow[i] == this->coocol[i] ) continue;
            coorow[inewpos] = this->coorow[i];
            coocol[inewpos] = this->coocol[i];
            cl1[inewpos] = this->cl1[i];
            cl2[inewpos] = this->cl2[i];
            fahl[inewpos] = this->fahl[i];
            fldr[inewpos] = this->fldr[i];
            //update the indptr array;
            int n = coorow[inewpos];
            if( n != lastn )
            {
                this->indptr[n] = inewpos;
                lastn = n;
            }
            inewpos += 1;
        }

        this->indptr[this->nodes] = inewpos;
        this->ncon = ncon;
        this->coorow = coorow;
        this->coocol = coocol;
        this->cl1 = cl1;
        this->cl2 = cl2;
        this->fahl = fahl;
        this->fldr = fldr;
    }

	/// data for Coo format 
	int ncon;
	int * coorow;
	int * coocol;
	int * indptr;
        
	/// total # of cells
	int nodes;
	
    /// size of nja storing distance to the shared of each connection
	float * cl1;
	
	/// size of nja storing reversed distance to the shared from each connection
	float * cl2;
	
	/// area of the connection of size nja
    float * fahl;
    
    /// direction indicator of size nja (+-1 for x, +-2 for y +-3 for z)
    int  * fldr;

    /// surface area from a plane view
    double * areas;
};


inline
ostream&  operator<<(ostream& out, CooData& data)
{
    out<<"nodes="<<data.nodes<<endl;
    out<<"ncon="<<data.ncon<<endl;
    out<<"coorow=(";
    for(int i=0;i<data.ncon;i++) out<<data.coorow[i]<<",";
    out<<")"<<endl;

    out<<"coocol=(";
    for(int i=0;i<data.ncon;i++) out<<data.coocol[i]<<",";
    out<<")"<<endl;

    out<<"area=(";
    for(int i=0;i<data.nodes;i++) out<<data.areas[i]<<",";
    out<<")"<<endl;

    out<<"cl1=(";
    for(int i=0;i<data.ncon;i++) out<<data.cl1[i]<<",";
    out<<")"<<endl;

    out<<"cl2=(";
    for(int i=0;i<data.ncon;i++) out<<data.cl2[i]<<",";
    out<<")"<<endl;

    out<<"fahl=(";
    for(int i=0;i<data.ncon;i++) out<<data.fahl[i]<<",";
    out<<")"<<endl;

    out<<"fldr=(";
    for(int i=0;i<data.ncon;i++) out<<data.fldr[i]<<",";
    out<<")"<<endl;

    return out;
}

} //namespace cusg
