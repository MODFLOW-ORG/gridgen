


#pragma once

#include "def_struct.h"
#include "CsrData.h"

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
    DefinitionExporter(){ m_number_per_line=10; }
    virtual ~DefinitionExporter(){}

    //export data to a definition file
    template<typename T>
	bool save(const string& filename, T& data)
	{
		ofstream fout(filename.c_str());
		if(fout.good()==false)
		{
			cerr<<"! Error: Failed to open file: "<<filename<<endl;
			return false; 
		}
		
		save(fout,data);
		fout.close();
		cout<<"- Save data to file: "<<filename<<endl;
		return true;
	}

    //

    template<typename T>
    bool save(const string& filename, T& data, int size)
    {
        ofstream fout(filename.c_str());
        if(fout.good()==false)
        {
            cerr<<"! Error: Failed to open file: "<<filename<<endl;
            return false;
        }
        save(fout,data,size);
        fout.close();
        cout<<"- Save data to file: "<<filename<<endl;
        return true;
    }
	
	//
    template<typename T>
    bool save(const string& filename, T& data, int size, int * formating)
    {
        ofstream fout(filename.c_str());
        if(fout.good()==false)
        {
            cerr<<"! Error: Failed to open file: "<<filename<<endl;
            return false;
        }
        save(fout,data,size,formating);
        fout.close();
        cout<<"- Save data to file: "<<filename<<endl;
        return true;
    }
    //

	bool save(ostream& out, ModflowGrid * grid);
	bool save(ostream& out, QuadTree3D * grid);
	bool save(ostream& out, grid_raw_data& data);
	bool save(ostream& out, modflow_grid_raw_data& data);
    bool save(ostream& out, quadtree_builder_raw_data& data);

    //
    bool save(grid_to_usgdata_raw_data& data);
    bool save(grid_intersection_raw_data& data);

    //save csr data into a set of files with "prefix"
    bool save(const string& prefix, CSRData& csr, int *nperlay);

    //setting flags
    void set_record_per_line(int i){
        if(i<1) i=1;
        if(i>2000) i=2000; //this is supposed to be FORTRAN limit
        m_number_per_line=i;
    }

    //save the grid
    bool save_node_coordinate(const string& prefix, Grid * grid);

protected:

    template<typename T>
    void save(ostream& out, const vector<T>& A)
    {
        int size=A.size();
        for(int i=0;i<size;i++)
        {
            save(out,A[i]);
            out<<" ";
        }
    }

    template<typename T>
    void save(ostream& out, const list<T>& A)
    {
        typedef typename list<T>::iterator IT;
        for(IT i=A.begin();i!=A.end();i++)
        {
            T& t=*i;
            save(out,t);
            out<<" ";
        }
    }

    //save an array of things
    template<typename T>
    void save(ostream& out, T * A, int size)
    {
        out.precision(10); //find a way to preserve the digit, adaptively
        for(int i=0;i<size;i++){
            if(i%m_number_per_line==0 && i>0) out<<"\n";
            out<<A[i]<<" ";
        }
    }

    //save an array of things with formating info
    template<typename T>
    void save(ostream& out, T * A, int size, int * formating)
    {
        int fid=0;
        int last_newline_id=0;
        out.precision(10); //find a way to preserve the digit, adaptively
        for(int i=0;i<size;i++){

            //check if a new formating line is reached
            if(i==formating[fid]-1){
                fid++;
                if(i>0 && i!=size-1){ last_newline_id=i; out<<"\n"; }
            }

            //the line is too long
            if( (i-last_newline_id)%m_number_per_line==0 && (i-last_newline_id)!=0 ) out<<"\n";

            //output the value
            out<<A[i]<<" ";
        }
    }

    void save(ostream& out, const double& data){ out<<data; }
    void save(ostream& out, const float& data){ out<<data; }
    void save(ostream& out, const int& data){ out<<data; }
    void save(ostream& out, const double * data){ out<<*data; }
    void save(ostream& out, float * data){ out<<*data; }
    void save(ostream& out, int * data){ out<<*data; }

private:
    int m_number_per_line; //
};

}
