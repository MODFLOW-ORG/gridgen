
/*! \file CsrData.h

\brief A structure to hold a sparse matrix in CSR format

\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/


#pragma once

#include <iostream>
using namespace std;

namespace cusg
{

/*! \brief Sparse Matrix in CSR format

    This object contains all of the unstructured grid information
    necessary for an unstructured grid model.  The data are stored
    in the CSR format.
 */
struct CSRData
{
	CSRData()
	{
	    nodes=0;
        nja=0;
        ja=NULL;
        ia=NULL;
        area=NULL;
        cl1=NULL;
        cl2=NULL;
        fahl=NULL;
        fldr=NULL;
        iac=NULL;
	}

	CSRData(int _nodes, int _nja, double * _area, int * _ia, int * _ja, float * _cl1, float * _cl2, float * _fahl, int * _fldr, int * _iac=NULL)
    {
        nodes=_nodes;
        nja=_nja;
        ja=_ja;
        ia=_ia;
        area=_area;
        cl1=_cl1;
        cl2=_cl2;
        fahl=_fahl;
        fldr=_fldr;
        iac=_iac;
    }

	void destory()
	{
	    delete [] ja;
	    delete [] ia;
	    delete [] area;
	    delete [] cl1;
	    delete [] cl2;
	    delete [] fahl;
	    delete [] fldr;
	    delete [] iac;
	}

	/// total # of cells
	int nodes;
	
	/// total # of connections
	int nja;
	
	/// array of size nja for connections
	int * ja;
	
	/// index array of size  nodes+1
	int * ia;
	
	/// surface area from a plane view, size is nodes
    double * area;
    
    /// size of nja storing distance to the shared face of each connection
	float * cl1;
	
	/// size of nja storing reversed distance to the shared face from each connection
	float * cl2;
	
	/// area of the connection of size nja
    float * fahl;
    
    /// direction indicator of size nja (+-1 for x, +-2 for y +-3 for z)
    int  * fldr;

    /// number of connection +1 for each cell, array of size  nodes
    int * iac;
};

inline
ostream&  operator<<(ostream& out, CSRData& data)
{
    out<<"nodes="<<data.nodes<<endl;
    out<<"nja="<<data.nja<<endl;
    out<<"ja=(";
    for(int i=0;i<data.nja;i++) out<<data.ja[i]<<",";
    out<<")"<<endl;

    out<<"ia=(";
    for(int i=0;i<data.nodes+1;i++) out<<data.ia[i]<<",";
    out<<")"<<endl;

    out<<"area=(";
    for(int i=0;i<data.nodes;i++) out<<data.area[i]<<",";
    out<<")"<<endl;

    out<<"cl1=(";
    for(int i=0;i<data.nja;i++) out<<data.cl1[i]<<",";
    out<<")"<<endl;

    out<<"cl2=(";
    for(int i=0;i<data.nja;i++) out<<data.cl2[i]<<",";
    out<<")"<<endl;

    out<<"fahl=(";
    for(int i=0;i<data.nja;i++) out<<data.fahl[i]<<",";
    out<<")"<<endl;

    out<<"fldr=(";
    for(int i=0;i<data.nja;i++) out<<data.fldr[i]<<",";
    out<<")"<<endl;

    return out;
}

inline
istream&  operator>>(istream& ins, CSRData& data)
{
    ins>>data.nodes;
    ins>>data.nja;

    data.ia = new int[data.nodes+1];
    data.ja = new int[data.nja];
    data.area = new double[data.nodes];
    data.cl1 = new float[data.nja];
    data.cl2 = new float[data.nja];
    data.fahl = new float[data.nja];
    data.fldr = new int[data.nja];

    assert(data.ia && data.ja && data.cl1 && data.cl2 && data.fahl && data.fldr);


    //locate space for CSR data here
    for(int i=0;i<data.nja;i++) ins>>data.ja[i];
    for(int i=0;i<data.nodes+1;i++) ins>>data.ia[i];
    for(int i=0;i<data.nodes;i++) ins>>data.area[i];
    for(int i=0;i<data.nja;i++) ins>>data.cl1[i];
    for(int i=0;i<data.nja;i++)  ins>>data.cl2[i];
    for(int i=0;i<data.nja;i++) ins>>data.fahl[i];
    for(int i=0;i<data.nja;i++) ins>>data.fldr[i];

    return ins;
}

} //namespace cusg
