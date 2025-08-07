#include "def_parser.h"

namespace cusg
{


//template<typename T> block_raw_data* createBlockType(block_raw_data& block) { return new T(block); }

//-----------------------------------------------------------------------------
//
//
// DefinitionParser
//
//
//-----------------------------------------------------------------------------

//constructor
DefinitionParser::DefinitionParser()
{
    m_block_type_map["modflow_grid"]     = &createBlockType<modflow_grid_raw_data>;
    m_block_type_map["quadtree"] = &createBlockType<quadtree_grid_raw_data>;
    m_block_type_map["quadtree_builder"] = &createBlockType<quadtree_builder_raw_data>;
    m_block_type_map["refinement_features"] = &createBlockType<refinement_raw_data>;
    m_block_type_map["active_domain"] = &createBlockType<active_domain_raw_data>;
    m_block_type_map["grid_to_shapefile"] = &createBlockType<gridtoshape_raw_data>;
	m_block_type_map["grid_to_usgdata"] = &createBlockType<grid_to_usgdata_raw_data>;
	m_block_type_map["grid_intersection"] = &createBlockType<grid_intersection_raw_data>;
	m_block_type_map["grid_to_vtkfile"] = &createBlockType<gridtovtk_raw_data>;
}

//desctructor
DefinitionParser::~DefinitionParser()
{
    typedef map<string, block_raw_data *>::iterator IT;
    for(IT i=m_blocks.begin();i!=m_blocks.end();i++)
    {
        delete i->second;
    }

    m_blocks.clear();
    m_block_type_map.clear();
}

bool DefinitionParser::parse(string filename)
{
    ifstream fin(filename.c_str());
    if(!fin.good()){
        cerr<<"! Error: Cannot open file:"<<filename<<endl;
        return false;
    }//end good

    bool r=parse(fin);
    fin.close();
    return r;
}

//parsing the ws file
bool DefinitionParser::parse(istream& in, string folder)
{
    while(!in.eof())
    {
        //parse raw data
        block_raw_data block;
        bool done=readBlock(in,block);
        if(done==false) break;

		block.folderPath = folder;//whether it is a relative path

        if(m_blocks.find(block.name)!=m_blocks.end()) //name exists
        {
            cerr<<"! Warning: duplicated record "<<block.name
                <<"found; replace the previous record"<<endl;
            delete m_blocks[block.name];
            m_blocks[block.name]=NULL;
        }

		//load the related path's .dfn file
		if(block.type == "load")
		{
			 string dfnPath = folder;
			 dfnPath += "/";
			 if(folder.empty() || folder =="")
				 dfnPath = "";

			 dfnPath += block.name;

			 //cout<<"dfnPath: "<<dfnPath<<endl;
			 int pathsize = dfnPath.size();
			 for (int i = 0; i < pathsize; i++)
			 {
				 if(dfnPath[i] =='\\')
					 dfnPath[i] = '/';
			 }
			 

			 ifstream relfin(dfnPath.c_str());
			 if(!relfin.good()){
				 cerr<<"! Error: Cannot open file:"<<dfnPath<<endl;
				 continue;
			 }//end good()

			  int dotIdx = dfnPath.find_last_of('/');
			  string relFolder ="";
			  if(dotIdx!=-1)
				  relFolder= dfnPath.substr(0, dotIdx);
			  
			  //cout<<"refFolder:"<<endl;
			  
			  bool r = parse(relfin, relFolder);
			  relfin.close();

			  continue;
		}

        //create block data with specific type
        if(m_block_type_map.find(block.type)!=m_block_type_map.end())
        {
            block_raw_data * ptr_block=m_block_type_map[block.type](block);
            if(ptr_block!=NULL)
            {
                bool r=ptr_block->parse();
                if(r) m_blocks[block.name]=ptr_block;
                else{
                    cerr<<"! Error: Something wrong paring the definition file"<<endl;
                    delete ptr_block;
                }
            }//end if
        }
        else
        {

            cerr<<"! Unknown block type:"<<block.type<<" ignored;\n"
                <<"! Do you forget to register block in DefinitionParser::DefinitionParser() "<<endl;

            continue;
        }//end if

    }//end while

    return true;
}

}//end of namespace


