

#pragma once

#include "def_struct.h"
#include "CsrData.h"
#include "RoundTrip.h"

namespace cusg
{

class Grid; //defined in Grid.h
class ModflowGrid;
class QuadTree3D;

//-----------------------------------------------------------------------------
//This class exports definition files

class DefinitionExporter
{
public:
    //constructor
    DefinitionExporter() { m_number_per_line = 10; }
    virtual ~DefinitionExporter() {}

    //export data to a definition file
    template <typename T> bool save(const string & filename, T & data)
    {
        ofstream fout(filename.c_str());
        if (fout.good() == false) throw GridgenError("cannot write file " + filename);

        save(fout, data);
        fout.close();
        cout << "- Save data to file: " << filename << endl;
        return true;
    }

    //

    template <typename T> bool save(const string & filename, const T & data, int size)
    {
        ofstream fout(filename.c_str());
        if (fout.good() == false) throw GridgenError("cannot write file " + filename);
        save(fout, data, size);
        fout.close();
        cout << "- Save data to file: " << filename << endl;
        return true;
    }

    //
    template <typename T>
    bool save(const string & filename, const T & data, int size, int * line_ends)
    {
        ofstream fout(filename.c_str());
        if (fout.good() == false) throw GridgenError("cannot write file " + filename);
        save(fout, data, size, line_ends);
        fout.close();
        cout << "- Save data to file: " << filename << endl;
        return true;
    }
    //

    bool save(ostream & out, QuadTree3D * grid);
    bool save(ostream & out, ModflowGridBlock & data);
    bool save(ostream & out, QuadtreeBuilderBlock & data);

    //
    bool save(GridToUsgdataBlock & data);

    //save csr data into a set of files with "prefix"; the face angles are for a grid rotated
    //by rotation_degrees
    bool save(const string & prefix, CSRData & csr, int * nperlay, double rotation_degrees);

    //setting flags
    void set_record_per_line(int i)
    {
        if (i < 1) i = 1;
        if (i > 2000) i = 2000; //this is supposed to be FORTRAN limit
        m_number_per_line = i;
    }

    //save the grid
    bool save_node_coordinate(const string & prefix, Grid * grid);

protected:
    template <typename T> void save(ostream & out, const vector<T> & a)
    {
        int size = a.size();
        for (int i = 0; i < size; i++)
        {
            save(out, a[i]);
            out << " ";
        }
    }

    template <typename T> void save(ostream & out, const list<T> & a)
    {
        using IT = typename list<T>::iterator;
        for (IT i = a.begin(); i != a.end(); i++)
        {
            T & t = *i;
            save(out, t);
            out << " ";
        }
    }

    //save an array of things
    template <typename T> void save(ostream & out, T * a, int size)
    {
        for (int i = 0; i < size; i++)
        {
            if (i % m_number_per_line == 0 && i > 0) out << "\n";
            out << round_trip_value(a[i]) << " ";
        }
    }

    //save an array of things with formatting info
    template <typename T> void save(ostream & out, T * a, int size, int * line_ends)
    {
        int fid = 0;
        int last_newline_id = 0;
        for (int i = 0; i < size; i++)
        {

            //check if a new formatting line is reached
            if (i == line_ends[fid] - 1)
            {
                fid++;
                if (i > 0 && i != size - 1)
                {
                    last_newline_id = i;
                    out << "\n";
                }
            }

            //the line is too long
            if ((i - last_newline_id) % m_number_per_line == 0 && (i - last_newline_id) != 0)
                out << "\n";

            //output the value
            out << round_trip_value(a[i]) << " ";
        }
    }

    void save(ostream & out, const double & data) { out << RoundTrip(data); }
    void save(ostream & out, const float & data) { out << RoundTrip(data); }
    void save(ostream & out, const int & data) { out << data; }
    void save(ostream & out, const double * data) { out << RoundTrip(*data); }
    void save(ostream & out, float * data) { out << RoundTrip(*data); }
    void save(ostream & out, int * data) { out << *data; }

private:
    int m_number_per_line; //
};

}
