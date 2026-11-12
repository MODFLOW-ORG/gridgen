
/*! \file CsrData.h

\brief A structure to hold a sparse matrix in CSR format

\author <a href="http://masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
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
    CSRData() : nodes(0), nja(0) {}

    CSRData(int nodes, int nja, vector<double> area, vector<int> ia, vector<int> ja,
            vector<double> cl1, vector<double> cl2, vector<double> fahl, vector<double> hwva,
            vector<int> fldr, vector<int> iac = vector<int>())
        : nodes(nodes), nja(nja), ja(std::move(ja)), ia(std::move(ia)), area(std::move(area)),
          cl1(std::move(cl1)), cl2(std::move(cl2)), fahl(std::move(fahl)), hwva(std::move(hwva)),
          fldr(std::move(fldr)), iac(std::move(iac))
    {
    }

    /// total # of cells
    int nodes;

    /// total # of connections
    int nja;

    /// array of size nja for connections
    vector<int> ja;

    /// index array of size  nodes+1
    vector<int> ia;

    /// surface area from a plane view, size is nodes
    vector<double> area;

    /// size of nja storing distance to the shared face of each connection
    vector<double> cl1;

    /// size of nja storing reversed distance to the shared face from each connection
    vector<double> cl2;

    /// area of the connection of size nja
    vector<double> fahl;

    /// face width of a horizontal connection and face area of a vertical one, size nja
    /// (MODFLOW 6 HWVA)
    vector<double> hwva;

    /// direction indicator of size nja (+-1 for x, +-2 for y +-3 for z)
    vector<int> fldr;

    /// number of connection +1 for each cell, array of size  nodes; empty if not computed
    vector<int> iac;

    /// 1 for a horizontal connection and 0 otherwise, size nja (MODFLOW 6 IHC)
    vector<int> ihc() const
    {
        vector<int> values(nja, 0);
        for (int i = 0; i < nja; i++)
            if (fldr[i] != 0 && std::abs(fldr[i]) != 3) values[i] = 1;
        return values;
    }

    /// 1 for a vertical connection and 0 otherwise, size nja (MODFLOW-USG IVC)
    vector<int> ivc() const
    {
        vector<int> values(nja, 0);
        for (int i = 0; i < nja; i++)
            if (std::abs(fldr[i]) == 3) values[i] = 1;
        return values;
    }

    /// angle of the outward normal of each horizontal face from the x axis, in degrees from 0
    /// to 360, for a grid rotated by rotation_degrees; 0 otherwise, size nja (MODFLOW 6 ANGLDEGX)
    vector<double> angldegx(double rotation_degrees) const
    {
        vector<double> values(nja, 0.0);
        for (int i = 0; i < nja; i++)
        {
            //+1 is +x, +2 is +y, -1 is -x, and -2 is -y in the grid's own axes
            double angle;
            if (fldr[i] == 1)
                angle = 0.0;
            else if (fldr[i] == 2)
                angle = 90.0;
            else if (fldr[i] == -1)
                angle = 180.0;
            else if (fldr[i] == -2)
                angle = 270.0;
            else
                continue;
            angle = std::fmod(angle + rotation_degrees, 360.0);
            if (angle < 0.0) angle += 360.0;
            //a negative angle closer to 0 than the rounding of 360 becomes 360
            if (angle >= 360.0) angle = 0.0;
            values[i] = angle;
        }
        return values;
    }
};

inline ostream & operator<<(ostream & out, CSRData & data)
{
    out << "nodes=" << data.nodes << endl;
    out << "nja=" << data.nja << endl;
    out << "ja=(";
    for (int i = 0; i < data.nja; i++) out << data.ja[i] << ",";
    out << ")" << endl;

    out << "ia=(";
    for (int i = 0; i < data.nodes + 1; i++) out << data.ia[i] << ",";
    out << ")" << endl;

    out << "area=(";
    for (int i = 0; i < data.nodes; i++) out << data.area[i] << ",";
    out << ")" << endl;

    out << "cl1=(";
    for (int i = 0; i < data.nja; i++) out << data.cl1[i] << ",";
    out << ")" << endl;

    out << "cl2=(";
    for (int i = 0; i < data.nja; i++) out << data.cl2[i] << ",";
    out << ")" << endl;

    out << "fahl=(";
    for (int i = 0; i < data.nja; i++) out << data.fahl[i] << ",";
    out << ")" << endl;

    out << "hwva=(";
    for (int i = 0; i < data.nja; i++) out << data.hwva[i] << ",";
    out << ")" << endl;

    out << "fldr=(";
    for (int i = 0; i < data.nja; i++) out << data.fldr[i] << ",";
    out << ")" << endl;

    return out;
}

inline istream & operator>>(istream & ins, CSRData & data)
{
    ins >> data.nodes;
    ins >> data.nja;

    data.ia.resize(data.nodes + 1);
    data.ja.resize(data.nja);
    data.area.resize(data.nodes);
    data.cl1.resize(data.nja);
    data.cl2.resize(data.nja);
    data.fahl.resize(data.nja);
    data.hwva.resize(data.nja);
    data.fldr.resize(data.nja);

    //locate space for CSR data here
    for (int i = 0; i < data.nja; i++) ins >> data.ja[i];
    for (int i = 0; i < data.nodes + 1; i++) ins >> data.ia[i];
    for (int i = 0; i < data.nodes; i++) ins >> data.area[i];
    for (int i = 0; i < data.nja; i++) ins >> data.cl1[i];
    for (int i = 0; i < data.nja; i++) ins >> data.cl2[i];
    for (int i = 0; i < data.nja; i++) ins >> data.fahl[i];
    for (int i = 0; i < data.nja; i++) ins >> data.hwva[i];
    for (int i = 0; i < data.nja; i++) ins >> data.fldr[i];

    return ins;
}

} //namespace cusg
