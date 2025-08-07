#include <ctype.h>
#include "def_struct.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include "QuadTreeBuilder.h"
#include "reader_treestruct.h"
#include "GridToShp.h"

namespace cusg
{

/* ********************************************************************** */
// skip blanks, tabs, line breaks and comment lines,
//  leaving us at the beginning of a token (or EOF)
//  (This code is taken from CORE lib)
inline int skip_comment_line (std::ifstream & in)
{
      int c;

      do {
        c = in.get();
        while ( c == '#' ) {
          do {// ignore the rest of this line
            c = in.get();
          } while ( c != '\n' );
          c = in.get(); // now, reach the beginning of the next line
        }
      } while (c == ' ' || c == '\t' || c == '\n' || c == '\r'); //ignore white spaces and newlines

      if (c == EOF)
        std::cout << "! Warning: Unexpected end of file." << std::endl;

      in.putback(c);  // this is non-white and non-comment char!
      return c;

}//skip_comment_line

// skips '\' followed by '\n'
//  NOTE: this assumes a very special file format (e.g., our BigInt File format)
//  in which the only legitimate appearance of '\' is when it is followed
//  by '\n' immediately!
//  (This code is taken from CORE lib)
inline int skip_backslash_new_line (std::istream & in)
{
      int c = in.get();

      while (c == '\\') {
        c = in.get();

        if (c == '\n' || c == '\r')
          c = in.get();
        else // assuming the very special file format noted above!
          cout<< "! Warning: continuation line \\ must be immediately followed by new line.\n";
      }//while
      return c;
}//skip_backslash_new_line


//-----------------------------------------------------------------------------
//
//
// BaseParser
//
//
//-----------------------------------------------------------------------------

bool BaseParser::readBlock(istream& in, block_raw_data& block)
{
    const int size=1024;
    const char * ignore=" \t[]()<>,!{}\r\n";
    char * tmp=new char[size];
    
    //
    while(!in.eof())
    {
        in.getline(tmp,size);

        if(strlen(tmp)==0) continue;
        //cout<<"get line="<<tmp<<endl;

        pair<char*,char*> name_value=split(tmp);

        list<string> name=tokenize(name_value.first,ignore);
        
        if( name.empty() )
        {
        	continue;
        }
        
        string label=name.front();
        if(label[0]=='#') continue; //comment

        label=tolower(label);

		/**********************************************************************/
		//here we simply regard "load" as a type 
		if(label == "load")
		{
			block.type="load";
			block.name = name.back();
			return true;
		}

        if(label=="begin")
        {
            name.pop_front();
            block.type=tolower(name.front());
            block.name=name.back();
        }
        else if(label=="end")
        {
            name.pop_front();
            string end_type=tolower(name.front());
        	if( block.type!=end_type )
        	{
        		cerr<<"! Error: block starts with type="<<block.type
        		    <<" but ends with a different type="<<end_type<<endl;
        	}
        	
        	break; //we are done for this block
        }
        else
        {
        	entry_raw_data data;
        	data.name=tolower(merge(name," "));

        	//cout<<"data.name="<<data.name<<endl;

			data.value=tokenize(name_value.second,ignore);
			block.push_back(data);
        }
    }//end while

    delete [] tmp;

	if(block.name.empty()) return false; //read nothing...

	return true;
}


//convert a string to a list of tokens
list<string> BaseParser::tokenize(char * tmp, const char * ignore)
{
    list<string> tokens;
    char * tok=strtok (tmp,ignore);
    while(tok!=NULL)
    {
        tokens.push_back(tok);
        tok=strtok(NULL,ignore);
    }
    return tokens;
}

// split a string into name and value by "="
pair<char*,char*> BaseParser::split(char * tmp)
{
	 char * pos=strchr(tmp,'=');
	 if(pos!=NULL){
	     *pos='\0';
	     pos++;
	     //cout<<"tmp="<<tmp<<", pos="<<pos<<"#"<<endl;
	 }

	 return make_pair(tmp,pos);
}

string BaseParser::merge(const list<string>& name, string del)
{
    string merged;
    typedef list<string>::const_iterator IT;
    IT last=--name.end();
    for(IT i=name.begin();i!=name.end();i++)
    {
        merged+=(*i);
        if(i!=last) merged+=del;
    }

    return merged;
}

string BaseParser::tolower(const string& s)
{
    string lower;
    for(string::const_iterator i=s.begin();i!=s.end();i++)
    {
        char c=::tolower(*i);
        lower.push_back(c);
    }
    return lower;
}
//-----------------------------------------------------------------------------
// split a string into name and value by "="
bool BaseParser::equal(const list<string>& name, string tmp)
{
	const char * ignore=" \t[]()<>,!{}\r\n";
	char * tmp2=new char[tmp.size()+1];
	assert(tmp2);
	strcpy(tmp2,tmp.c_str());
	list<string> name2=tokenize(tmp2, ignore);
	delete [] tmp2;
	if( name.size()!=name2.size() )  return false;
	//
	typedef list<string>::const_iterator IT;
	for(IT i=name.begin(), j=name2.begin();i!=name.end();i++,j++)
	{
		if( (*i)!=(*j) ) return false;
	}
	//
	return true;
}

//-----------------------------------------------------------------------------
//
//
// modflow_grid_raw_data
//
//
//-----------------------------------------------------------------------------

bool modflow_grid_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["rotation_angle"];
    if(!value.empty()) 
		{
			rotation_angle=atof(value.front().c_str());
			rotation_angle = rotation_angle * 2 *3.1415926 /360;
	}
    else{
        cerr<<"! Error: missing rotation_angle"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["x_offset"];
    if(!value.empty()) x_offset=atof(value.front().c_str());
    else{
        cerr<<"! Error: missing x_offset"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["y_offset"];
    if(!value.empty()) y_offset=atof(value.front().c_str());
    else{
        cerr<<"! Error: missing y_offset"<<endl;
        assert(false);
		return false;
    }

	
    value=name_value_map["length_unit"];
    if(!value.empty()) length_unit=value.front();
    else
	{
        //cerr<<"! Error: missing length_unit"<<endl;
        //assert(false);
		//return false;
    }

    value=name_value_map["nlay"];
    if(!value.empty()) nlay=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing nlay"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["nrow"];
    if(!value.empty()) nrow=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing nrow"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["ncol"];
    if(!value.empty()) ncol=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing ncol"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["delr"];
    if(!value.empty()){
        delr_str=value;
    }
    else{
        cerr<<"! Error: missing delr"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["delc"];
    if(!value.empty()){
        delc_str=value;
    }
    else{
        cerr<<"! Error: missing delc"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["top"];
    if(!value.empty()){
        top_str=value;
    }
    else{
        cerr<<"! Error: missing top layer"<<endl;
        assert(false);
		return false;
    }


    for(int i=0;i<nlay;i++)
    {
        char layername[36];
        sprintf(layername,"bottom layer %d",i+1);
        value=name_value_map[layername];
        if(!value.empty()){
            bot_str[i]=value;
			//relative path
        }
        else{
            cerr<<"! Error: missing bottom layer "<<i + 1 << "("<<i<<" for 0 based numbering)"<<endl;
            assert(false);
			return false;
        }
    }


    return true;
}

bool modflow_grid_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
    if(m_initizlied) return true;

    delr=createArray<REAL>(delr_str,ncol, getFolderPath());
    delc=createArray<REAL>(delc_str,nrow, getFolderPath());
    top=createArray<REAL>(top_str,ncol*nrow, getFolderPath());

    bot.reserve(nlay+1);
    bot.push_back(top);
    for(int i=0;i<nlay;i++)
    {
        REAL * B=createArray<REAL>(bot_str[i],ncol*nrow, getFolderPath());
        bot.push_back(B);
    }

    m_initizlied=true;

    return true;
}

Grid * modflow_grid_raw_data::build()
{
    if(mfgrid==NULL)
    {
        mfgrid=new ModflowGrid(*this);
        assert(mfgrid);
        mfgrid->name=this->name;
    }

    return (Grid *)mfgrid;
}

//-----------------------------------------------------------------------------
//
//
// quadtree_grid_raw_data
//
//
//-----------------------------------------------------------------------------

bool quadtree_grid_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["modflow_grid"];
    if(!value.empty()) modflow_grid=value.front();
    else{
        cerr<<"! Error: missing modflow_grid"<<endl;
        assert(false);
		return false;
    }

	value = name_value_map["structure_file"];
	if(!value.empty()) tree_structure_file = value.back();//value.front();
	else
	{
		cerr<<"! Error: missing tree_structure_file"<<endl;
		assert(false);
		return false;
	}
	
    for(int i=0;;i++)
    {
        char layername[36];
        //
        sprintf(layername,"top layer %d",i+1);
        value=name_value_map[layername];
        if(!value.empty()) top_str[i]=value;
        //
        sprintf(layername,"bottom layer %d",i+1);
        value=name_value_map[layername];
        if(!value.empty()) bot_str[i]=value;
        else break;
        //
    }

    return true;
}


bool quadtree_grid_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
    if(m_initizlied) return true;

    this->modflow_data=dynamic_cast<modflow_grid_raw_data *>(blocks[modflow_grid]);

    if(modflow_data==NULL)
    {
        cerr<<"! Error: Cannot resolve referenced modflow_grid name: "<<modflow_grid<<endl;
        return false;
    }

    if(modflow_data->isInitialized()==false)
    {
       if( modflow_data->initialize(blocks)==false )
           return false;
    }

    //read/open this->tree_structure_file here
	parseTreeStructFile(tree_structure_file,allNodes, getFolderPath());

//    nodesperlay=createArray<int>(nodesperlay_str,modflow_data->nlay, getFolderPath());
//    ia=createArray<int>(ia_str,this->nodes+1, getFolderPath());
//    ja=createArray<int>(ja_str,this->ncon, getFolderPath());
//    connection_direction_array=createArray<int>(connection_direction_array_str,this->ncon, getFolderPath());

	/*
    top.reserve(modflow_data->nlay);
    bot.reserve(modflow_data->nlay);
    
	//int ncol_nrow=modflow_data->ncol*modflow_data->nrow; //TODO: this is incorrect...
    
	for(int i=0;i<modflow_data->nlay;i++)
    {

        REAL * T=createArray<REAL>(top_str[i],ncol_nrow, getFolderPath());
        REAL * B=createArray<REAL>(bot_str[i],ncol_nrow, getFolderPath());
        top.push_back(T);
        bot.push_back(B);
    }
	*/

    m_initizlied=true;
    return true;
}

Grid * quadtree_grid_raw_data::build()
{
    //
    if(this->modflow_data==NULL)
    {
        cerr<<"! Error: Cannot resolve referenced modflow_grid name: "<<modflow_grid<<endl;
        return NULL;
    }

    if(qtree_grid==NULL) //if it's not already built
    {
        ModflowGrid * modflowGrid=dynamic_cast<ModflowGrid*>(modflow_data->build());
		if(modflowGrid==NULL)
		{
			cerr<<"! Error: quadtree_grid_raw_data::build: Cannot retrieve ModflowGrid"<<endl;
			return NULL; 
		}

        qtree_grid=new QuadTree3D(modflowGrid);
        assert(qtree_grid);

        //fill the QuadTree3d
        fillQuadTree3D();

		//copy surface data from top/bottom to those in qtree_grid
		int nodesize=qtree_grid->nodes();
		if(qtree_grid->top!=NULL) delete [] qtree_grid->top;
		if(qtree_grid->bot!=NULL) delete [] qtree_grid->bot;
		qtree_grid->top=new double[nodesize];
		qtree_grid->bot=new double[nodesize];
		assert(qtree_grid->top && qtree_grid->bot);
		int index=0;
		for(int layer=0;layer<modflowGrid->nlay;layer++)
		{
			int node_per_layer = qtree_grid->nodegroup.nodelay_array[layer];
			
			double * top_tmp=createArray<REAL>(top_str[layer],node_per_layer, getFolderPath());
			memcpy( &qtree_grid->top[index], top_tmp, sizeof(double)*node_per_layer );
			delete [] top_tmp;//free data

			double * bot_tmp=createArray<REAL>(bot_str[layer],node_per_layer, getFolderPath());
			memcpy( &qtree_grid->bot[index], bot_tmp, sizeof(double)*node_per_layer );
			delete [] bot_tmp;//free data

			index+=node_per_layer;
		}//end for layer

		//done
        qtree_grid->name=this->name;
    }

    return (Grid *)qtree_grid;
}

void quadtree_grid_raw_data::fillQuadTree3D()
{
	int count = 0;
	for (vector<TreeStructNode>::iterator nite = allNodes.begin(); nite != allNodes.end(); ++nite, count++)
	{

		TreeStructNode& tmpNode = *nite;
		int layerNum = tmpNode.layer - 1;
		int row = tmpNode.row - 1;
		int col = tmpNode.col - 1;

//		Index2d id2d=getIndex(qtree_grid->getModflowGrid(),test);
//		int id=mf->get_nodeid(layerNum,id2d[1],id2d[0]);
//		Box * box=grid->nodegroup[id];

		
		int nodeid = qtree_grid->getModflowGrid()->get_nodeid(layerNum, row, col);
		Box* nodeobj = qtree_grid->nodegroup[nodeid];

		if( (strcmp(tmpNode.code.c_str(),"")== 0) || tmpNode.code.empty())
		{
			nodeobj->number = tmpNode.number;
			nodeobj->isLeaf = true;
			if(nodeobj->number==-1) nodeobj->active=false;
		}
		else
		{
			Box* tmpNodeBoxObj = nodeobj;
			int codeDepth = tmpNode.code.size();
			int depthIdx = 0;
			while(depthIdx < codeDepth)
			{
				if(tmpNodeBoxObj->isLeaf)
				{
					bool r = tmpNodeBoxObj->split(0);
					assert(r);

					for (int c = 0;c < 4; c++)
					{
					    int id=c;
					    if(c==3) id=2;
					    if(c==2) id=3;
						tmpNodeBoxObj->pChildren[id]->isLeaf = true;
						tmpNodeBoxObj->pChildren[id]->updateStatus();
						qtree_grid->nodegroup.add(tmpNodeBoxObj->pChildren[id]);
					}
					tmpNodeBoxObj->isLeaf = false;
				}
				char t = tmpNode.code[depthIdx];
				switch (t)
				{
				case '1':
					{
						tmpNodeBoxObj = tmpNodeBoxObj->pChildren[0];
						break;
					}
				case '2':
					{
						tmpNodeBoxObj = tmpNodeBoxObj->pChildren[1];
						break;
					}
				case '3':
					{
					    tmpNodeBoxObj = tmpNodeBoxObj->pChildren[3];
						break;
					}
				case '4':
					{
					    tmpNodeBoxObj = tmpNodeBoxObj->pChildren[2];
						break;
					}
				default:
					{
						cerr<<"! Error: splitting box error!"<<endl;
						break;
					}
				}
				depthIdx++;
			}

			tmpNodeBoxObj->number = tmpNode.number;
			tmpNodeBoxObj->isLeaf = true;
			if(tmpNodeBoxObj->number==-1) tmpNodeBoxObj->active=false;

			//if(!(tmpNodeBoxObj->isLeaf && tmpNodeBoxObj->number!=0))
			//	cout<<"error"<<endl;
			//count++;
		}
	}
	//cout<<"total node number "<<count<<endl;

	//refresh the tree to get the new node id array
	qtree_grid->refresh();

	/********************************************************************/
	//check whether the quadtree is built successfully
#if DEBUG
	cout<<"********************************************************"<<endl;
	int size = qtree_grid->nodegroup.size();
	cout<<"The total number of size is : "<<size<<endl;
	int showNum = 0;
	for(int i=0;i<size;i++)
	{
		Box * box=qtree_grid->nodegroup[i];
		if(box->isLeaf){
			showNum++;
		}
	}
	cout<<"The leaf nodes' number is : "<<showNum<<endl;
	cout<<"********************************************************"<<endl;
#endif

}
//-----------------------------------------------------------------------------
//
//
// quadtree_builder_raw_data
//
//
//-----------------------------------------------------------------------------

bool quadtree_builder_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["max_smoothing_level_target"];
    if(!value.empty()) smoothing_level_vertical=smoothing_level_horizontal=atoi(value.front().c_str());

	value=name_value_map["smoothing_level_horizontal"];
    if(!value.empty()) smoothing_level_horizontal=atoi(value.front().c_str());

	value=name_value_map["smoothing_level_vertical"];
    if(!value.empty()) smoothing_level_vertical=atoi(value.front().c_str());

    value=name_value_map["modflow_grid"];
    if(!value.empty()) modflow_grid=value.front();
    else{
        cerr<<"! Error: missing modflow_grid"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["grid_definition_file"];
    if(!value.empty()) grid_definition_file=value.front();
    else{
        cerr<<"! Error: missing grid_definition_file"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["smoothing"];
    if(!value.empty()) smoothing=(tolower(value.front())=="full")?true:false;

    value=name_value_map["node_numbering"];
    if(!value.empty()){
        if(tolower(value.front())=="zero-based")
            one_based_node_numbering=false;
    }

	value=name_value_map["autoalignment"];
	if(!value.empty()) auto_alignment=(tolower(value.front())=="true")?true:false;

    //parsing active domain
    for(int i=0;;i++)
    {
        char activename[64];
        sprintf(activename,"active_domain layer %d",i+1);
        value=name_value_map[activename];
        if(!value.empty()) active_domain_str[i]=value;
        else break;
    }

    //parsing refinement stuff
    for(int i=0;;i++)
    {
        char layername[64];

        //get info about refinement features
        sprintf(layername,"refine_by_features layer %d",i+1);
        value=name_value_map[layername];
        if(!value.empty()) refine_by_features_str[i]=value;
        else //we also allow refinement_features
        {
			//
			sprintf(layername,"refinement_features layer %d",i+1);
			value=name_value_map[layername];
			if(!value.empty()) refine_by_features_str[i]=value;
        }
        
        //get info about refinement array
        sprintf(layername,"refine_by_array layer %d",i+1);
        value=name_value_map[layername];
        if(!value.empty()) refine_by_array_str[i]=value;
        else
		{
			sprintf(layername,"refinement_array layer %d",i+1);
			value=name_value_map[layername];
			if(!value.empty()) refine_by_array_str[i]=value;        
        }
        
        //surface array for each layer
        sprintf(layername,"top layer %d",i+1);
        value=name_value_map[layername];
        if(!value.empty()) top_str[i]=value;
        //
        sprintf(layername,"bottom layer %d",i+1);
        value=name_value_map[layername];

        //
        if(!value.empty()) bot_str[i]=value;
        else break;
        //
    }

    return true;
}

bool quadtree_builder_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
    if(m_initizlied) return true;

    //initialize modflow_data if it has not been initialized yet
    this->modflow_data=dynamic_cast<modflow_grid_raw_data *>(blocks[modflow_grid]);
    if(modflow_data==NULL)
    {
        cerr<<"! Error: InQuadtree Builder: Cannot resolve referenced modflow_grid name: "<<modflow_grid<<endl;
        return false;
    }

    if(modflow_data->isInitialized()==false)
    {
       if( modflow_data->initialize(blocks)==false )
           return false;
    }

    // initialize return false;
    active_domain_data.reserve(modflow_data->nlay);
    for(int i=0;i<modflow_data->nlay;i++)
    {
        list<active_domain_raw_data*> aF;
        list<string>& aF_str=active_domain_str[i];
        for(list<string>::iterator s=aF_str.begin();s!=aF_str.end();s++)
        {
            active_domain_raw_data* af=dynamic_cast<active_domain_raw_data *>(blocks[*s]);
            if(af!=NULL){
                if(af->initialize(blocks)==false)
                {
                    return false;
                }
                else aF.push_back(af);
            }
            else{
                cerr<<"! Error: Quadtree Builder: Block "<<*s<<" is referenced but not defined"<<endl;
				assert(false);
				return false;
            }
        }// end for s
        active_domain_data.push_back(aF);
    }

    // initialize refinement features
    refine_by_features.reserve(modflow_data->nlay);
    refine_by_array.reserve(modflow_data->nlay);
    top.reserve(modflow_data->nlay);
    bot.reserve(modflow_data->nlay);
    int ncol_nrow=modflow_data->ncol*modflow_data->nrow;

	//fill in missing info about the depth values of each layer
	for(int i=0;i<modflow_data->nlay;i++)
	{
		if(top_str[i].empty()){ top_str[i].push_back("replicate"); top_str[i].push_back(modflow_data->name); }
		if(bot_str[i].empty()){ bot_str[i].push_back("replicate"); bot_str[i].push_back(modflow_data->name); }
	}

	//process each layer now
    for(int i=0;i<modflow_data->nlay;i++)
    {
        //refine by features
        list<refinement_raw_data*> rF;
        list<string>& rF_str=refine_by_features_str[i];
        for(list<string>::iterator s=rF_str.begin();s!=rF_str.end();s++)
        {
            refinement_raw_data * rf=dynamic_cast<refinement_raw_data *>(blocks[*s]);
            if(rf!=NULL) rF.push_back(rf);
			else
			{
				cerr<<"! Error: Quadtree Builder: There is no refinement feature named "<<*s<<endl;
				assert(false);
				return false;
			}
        }

        refine_by_features.push_back(rF);

        //refine by array
        int *  R=createArray<int>(refine_by_array_str[i],ncol_nrow, getFolderPath());
        refine_by_array.push_back(R);

		//surface data
		top.push_back( make_pair( tolower(top_str[i].front()), top_str[i].back()) );
		bot.push_back( make_pair( tolower(bot_str[i].front()),bot_str[i].back()) );

		//interpolation mode
		if(top_str[i].size() == 3)
		{
			vector<string> vs(top_str[i].begin(), top_str[i].end());
			if(tolower(vs[1]) == "area_weighted")
				top_area_weighted.insert(make_pair(i, true));
			else if(tolower(vs[1]) == "no_area_weighted")
				top_area_weighted.insert(make_pair(i, false));
		}
		if(bot_str[i].size() == 3)
		{
			vector<string> vs(bot_str[i].begin(), bot_str[i].end());
			if(tolower(vs[1]) == "area_weighted")
				bot_area_weighted.insert(make_pair(i, true));
			else if(tolower(vs[1]) == "no_area_weighted")
				bot_area_weighted.insert(make_pair(i, false));
		}
    }

    m_initizlied=true;
    return true;
}

Grid * quadtree_builder_raw_data::build()
{
    //
    if(this->modflow_data==NULL)
    {
        cerr<<"! Error: Quadtree Builder: Cannot resolve referenced modflow_grid name: "<<modflow_grid<<endl;
        return NULL;
    }

    if(qtree_grid==NULL) //if it's not already built
    {
        Grid * modflowGrid=modflow_data->build();

        qtree_grid=new QuadTree3D((ModflowGrid*)modflowGrid);
        assert(qtree_grid);
        if(this->grid_definition_file.empty()) //no name defined...
            qtree_grid->name=this->name;
        else{
            int pos=this->grid_definition_file.find_last_of('.');
            qtree_grid->name=this->grid_definition_file.substr(0,pos);
        }

        Quadtree_Builder builder;
        builder.build(qtree_grid,*this);
    }

    return (Grid *)qtree_grid;
}

//-----------------------------------------------------------------------------
//
//
// surface_interpolation_raw_data
//
//
//-----------------------------------------------------------------------------

bool surface_interpolation_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["filename"];
    if(!value.empty()) filename=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing filename"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["surface_type"];
    if(!value.empty()) surface_type=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing surface_type"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["modflow_grid"];
    if(!value.empty()) modflow_grid=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing modflow_grid"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["interpolation_type"];
    if(!value.empty()) interpolation_type=atoi(value.front().c_str());
    else{
        cerr<<"! Error: missing interpolation_type"<<endl;
        assert(false);
		return false;
    }

    return true;
}

//-----------------------------------------------------------------------------
//
//
// grid_intersection_raw_data
//
//
//-----------------------------------------------------------------------------

bool grid_intersection_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["grid"];
    if(!value.empty()) grid_name=value.front();
    else{
        cerr<<"! Error: missing grid"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["shapefile"];
    if(!value.empty())
	{
			shapefile=value.front();
			//replace the '\\'
			string fullpath = getFolderPath();
			fullpath += shapefile;
			shapefile = fullpath;

			int pathsize = shapefile.size();
			for (int i = 0;i < pathsize; i++)
			{
				if(shapefile[i] =='\\')
					shapefile[i] = '/';
			}
	}
    else{
        cerr<<"! Error: missing shapefile"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["feature_type"];
    if(!value.empty()) feature_type=value.front();
    else{
        cerr<<"! Error: missing feature_type"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["output_file"];
    if(!value.empty())
	{
			output_file=value.front();
			string fullpath = getFolderPath();
			fullpath +=output_file;
			output_file = fullpath;
			
			int pathsize = output_file.size();
			for(int i = 0; i< pathsize; i++)
			{
				if(output_file[i] =='\\')
					output_file[i] = '/';
			}
	}
    else{
        cerr<<"! Error: missing output_file"<<endl;
        assert(false);
		return false;
    }

	value= name_value_map["output_shapefile"];
	if(!value.empty())
	{
			output_shapefile = value.front();
			string fullpath = getFolderPath();
			fullpath +=output_shapefile;
			
			output_shapefile = fullpath;
			int pathsize = output_shapefile.size();
			for (int i = 0; i < pathsize; i++)
			{
				if(output_shapefile[i] == '\\')
						output_shapefile[i] = '/';
			}
	}

	value = name_value_map["output_vtkfile"];
	if(!value.empty())
	{
			output_vtkfile = value.front();
			string fullpath = getFolderPath();
			fullpath +=output_vtkfile;
			
			output_vtkfile = fullpath;
			int pathsize = output_vtkfile.size();
			for (int i = 0; i < pathsize; i++)
			{
				if(output_vtkfile[i] == '\\')
						output_vtkfile[i] = '/';
			}
	}

    value=name_value_map["layer"];
    if(!value.empty()) layer=atoi(value.front().c_str())-1;
    else{
        cerr<<"! Error: missing layer"<<endl;
        assert(false);
		return false;
    }

	value = name_value_map["attributes"];
	if(!value.empty())
	{
		attributes.insert(attributes.end(), value.begin(), value.end());
	}
	else
	{
		attributes.push_back("all");
	}

    return true;
}

bool grid_intersection_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
    if(m_initizlied) return true;

    //get the grid
    block_raw_data * block=blocks[grid_name];

    grid_raw_data* grid=dynamic_cast<grid_raw_data*>(block);
    if( grid==NULL )
    {
        cerr<<"! Error: "<<grid_name<<" is not a grid"<<endl;
        return false;
    }

    //initialize the grid if necessary
    bool r=grid->initialize(blocks);
    if(r==false){
        cerr<<"! Error: initialize "<<grid_name<<" failed"<<endl;
        return false;
    }

    //should we check if the grid is built...
    //it's safe to built a grid multiple times anyway...
    this->grid=grid->build();

    //check again
    if(this->grid==NULL)
    {
        cerr<<"! Error: build "<<grid_name<<" failed"<<endl;
        return false;
    }

    m_initizlied=true;

    return true;
}

//-----------------------------------------------------------------------------
//
//
// active_domain_raw_data
//
//
//-----------------------------------------------------------------------------

bool active_domain_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
//    list<string>& value=name_value_map["grid"];
//    if(!value.empty()) grid_name=value.front();
//    else{
//        cerr<<"! Error: missing grid"<<endl;
//        assert(false);
//    }

    list<string>& value=name_value_map["shapefile"];
    if(!value.empty()) shapefile=value.front();
    else{
        cerr<<"! Error: missing shapefile"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["feature_type"];
    if(!value.empty()) feature_type=value.front();
    else{
        cerr<<"! Error: missing feature_type"<<endl;
        assert(false);
		return false;
    }

//    value=name_value_map["layer"];
//    if(!value.empty()){
//        for(list<string>::iterator i=value.begin();i!=value.end();i++)
//        {
//            layers.push_back( atoi(i->c_str()) );
//        }
//    }
//    else{
//        cerr<<"! Error: missing layer"<<endl;
//        assert(false);
//    }

    value=name_value_map["include_boundary"];
    if(!value.empty()) include_boundary=(::tolower(value.front()[0])=='f')?false:true;

    return true;
}

bool active_domain_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
    if(m_initizlied) return true;

//    //get the grid
//    block_raw_data * block=blocks[grid_name];
//
//    grid_raw_data* grid=dynamic_cast<grid_raw_data*>(block);
//    if( grid==NULL )
//    {
//        cerr<<"! Error: "<<grid_name<<" is not a grid"<<endl;
//        return false;
//    }
//
//    //initialize the grid if necessary
//    bool r=grid->initialize(blocks);
//    if(r==false){
//        cerr<<"! Error: initialize "<<grid_name<<" failed"<<endl;
//        return false;
//    }
//
//    //should we check if the grid is built...
//    //it's safe to built a grid multiple times anyway...
//    this->grid=grid->build();
//
//    //check again
//    if(this->grid==NULL)
//    {
//        cerr<<"! Error: build "<<grid_name<<" failed"<<endl;
//        return false;
//    }

    m_initizlied=true;

    return true;
}

//-----------------------------------------------------------------------------
//
//
// refinement_raw_data
//
//
//-----------------------------------------------------------------------------

bool refinement_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["shapefile"];
    if(!value.empty()) shapefile=value.front();
    else{
        cerr<<"! Error: missing shapefile"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["feature_type"];
    if(!value.empty()) featuretype=value.front();
    else{
        cerr<<"! Error: missing feature_type"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["refinement_level"];
    if(!value.empty()) refinement_level=atoi(value.front().c_str());
    else{
        value=name_value_map["refinement_level_by_attribute"];
        if(!value.empty()) refinement_level_by_attribute=value.front();
        else{
            cerr<<"! Missing refinement_level or refinement_level_by_attribute"<<endl;
            assert(false);
			return false;
        }
    }

//    value=name_value_map["refinement_level_by_attribute"];
//    if(!value.empty()) refinement_level_by_attribute=value.front();
//    else{
//        cerr<<"! Error: missing refinement_level_by_attribute"<<endl;
//        assert(false);
//    }

    return true;
}

//-----------------------------------------------------------------------------
//
//
// grid_to_usgdata_raw_data
//
//
//-----------------------------------------------------------------------------

bool grid_to_usgdata_raw_data::parse()
{
    map<string, list<string> > name_value_map;

    for(iterator i=begin(); i!=end(); i++)
    {
        string& name=i->name;
        list<string>& value=i->value;

        if( name.empty() || value.empty() )
        {
            cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
            continue;
        }

        name_value_map[name]=value;
    }


    //check the names
    list<string>& value=name_value_map["grid"];
    if(!value.empty()) grid_name=value.front();
    else{
        cerr<<"! Error: missing grid"<<endl;
        assert(false);
		return false;
    }

    value=name_value_map["usg_data_prefix"];
    if(!value.empty()) usg_data_prefix=value.front();
    else{
        cerr<<"! Error: missing usg_data_prefix"<<endl;
        assert(false);
		return false;
    }

	value = name_value_map["vertical_pass_through"];
	if(!value.empty())
	{
		string vs = tolower(value.front());
		if(vs == "true")
			vertical_pass_through = true;
	}

    return true;
}

bool grid_to_usgdata_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
    if(m_initizlied) return true;

    //get the grid
    block_raw_data * block=blocks[grid_name];

    grid_raw_data* grid=dynamic_cast<grid_raw_data*>(block);
    if( grid==NULL )
    {
        cerr<<"! Error: "<<grid_name<<" is not a grid"<<endl;
        return false;
    }

    //initialize the grid if necessary
    bool r=grid->initialize(blocks);
    if(r==false){
        cerr<<"! Error: initialize "<<grid_name<<" failed"<<endl;
        return false;
    }

    //should we check if the grid is built...
    //it's safe to built a grid multiple times anyway...
    this->grid=grid->build();

    //check again
    if(this->grid==NULL)
    {
        cerr<<"! Error: build "<<grid_name<<" failed"<<endl;
        return false;
    }

    m_initizlied=true;

    return true;
}

//-----------------------------------------------------------------------------
//
//
// gridtoshape_raw_data
//
//
//-----------------------------------------------------------------------------



bool gridtoshape_raw_data::parse()
{
	map<string, list<string> > name_value_map;

	for(iterator i=begin(); i!=end(); i++)
	{
        string& name=i->name;
		list<string>& value=i->value;

		if( name.empty() || value.empty() )
		{
			cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
			continue;
		}

		name_value_map[name]=value;
	}


	//check the names
	list<string>& value=name_value_map["grid"];
	if(!value.empty()) grid_name=value.front();
	else{
		cerr<<"! Error: missing grid"<<endl;
		assert(false);
		return false;
	}

	value = name_value_map["shapefile"];
	if(!value.empty())
	{		
			shapefile = value.back();//value.front();
			string fullpath = getFolderPath();
			fullpath += shapefile;
			shapefile = fullpath;

			int pathsize = shapefile.size();
			for(int i = 0; i < pathsize; i++)
			{
				if(shapefile[i] =='\\')
					shapefile[i] = '/';
			}
	}
	else
	{
		cerr<<"! Error: missing shapefile"<<endl;
		assert(false);
		return false;
	}

	value=name_value_map["node_numbering"];
    if(!value.empty()){
        if(tolower(value.front())=="zero-based")
            one_based_node_numbering=false;
    }

	value=name_value_map["without_inactive"];
    if(!value.empty()){
        if(tolower(value.front())=="false")
            without_inactive=false;
    }

	value=name_value_map["concide"];
    if(!value.empty()){
        if(tolower(value.front())=="false")
            concide=false;
    }

	value=name_value_map["feature_type"];
	//if(!value.empty()) feature_type=atoi(value.front().c_str());
	if(!value.empty()) feature_type = value.front();
	else{
		cerr<<"! Error: missing feature_type"<<endl;
		assert(false);
		return false;
	}

	return true;
}
bool gridtoshape_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
	if(m_initizlied) return true;

    //get the grid
    block_raw_data * block=blocks[grid_name];

    grid_raw_data* grid=dynamic_cast<grid_raw_data*>(block);
    if( grid==NULL )
    {
        cerr<<"! Error: "<<grid_name<<" is not a grid"<<endl;
        return false;
    }

    //initialize the grid if necessary
    bool r=grid->initialize(blocks);
    if(r==false){
        cerr<<"! Error: initialize "<<grid_name<<" failed"<<endl;
        return false;
    }

    //should we check if the grid is built...
    //it's safe to built a grid multiple times anyway...
    this->grid=grid->build();

    //check again
    if(this->grid==NULL)
    {
        cerr<<"! Error: build "<<grid_name<<" failed"<<endl;
        return false;
    }

	m_initizlied=true;
	return true;
}

//-----------------------------------------------------------------------------
//
//
// gridtovtk_raw_data
//
//
//-----------------------------------------------------------------------------


bool gridtovtk_raw_data::parse()
{
	map<string, list<string> > name_value_map;

	for(iterator i=begin(); i!=end(); i++)
	{
        string& name=i->name;
		list<string>& value=i->value;

		if( name.empty() || value.empty() )
		{
			cerr<<"! Error: Cannot parse record: empty name or empty value"<<endl;
			continue;
		}

		name_value_map[name]=value;
	}


	//check the names
	list<string>& value=name_value_map["grid"];
	if(!value.empty()) grid_name=value.front();
	else{
		cerr<<"! Error: missing grid"<<endl;
		assert(false);
		return false;
	}

	value = name_value_map["vtkfile"];
	if(!value.empty())
	{		
			vtkfile = value.back();//value.front();
			string fullpath = getFolderPath();
			fullpath += vtkfile;
			vtkfile = fullpath;

			int pathsize = vtkfile.size();
			for(int i = 0; i < pathsize; i++)
			{
				if(vtkfile[i] =='\\')
					vtkfile[i] = '/';
			}
	}
	else
	{
		cerr<<"! Error: missing vtkfile"<<endl;
		assert(false);
		return false;
	}
		
	value=name_value_map["without_inactive"];
    if(!value.empty()){
        if(tolower(value.front())=="false")
            without_inactive=false;
    }

	value=name_value_map["node_numbering"];
    if(!value.empty()){
        if(tolower(value.front())=="zero-based")
            one_based_node_numbering=false;
    }

	value=name_value_map["share_vertex"];
	if(!value.empty()) b_share_vertex=(tolower(value.front())=="true")?true:false;

	return true;
}
bool gridtovtk_raw_data::initialize(map<string, block_raw_data *>& blocks)
{
	if(m_initizlied) return true;

    //get the grid
    block_raw_data * block=blocks[grid_name];

    grid_raw_data* grid=dynamic_cast<grid_raw_data*>(block);
    if( grid==NULL )
    {
        cerr<<"! Error: "<<grid_name<<" is not a grid"<<endl;
        return false;
    }

    //initialize the grid if necessary
    bool r=grid->initialize(blocks);
    if(r==false){
        cerr<<"! Error: initialize "<<grid_name<<" failed"<<endl;
        return false;
    }

    //should we check if the grid is built...
    //it's safe to built a grid multiple times anyway...
    this->grid=grid->build();

    //check again
    if(this->grid==NULL)
    {
        cerr<<"! Error: build "<<grid_name<<" failed"<<endl;
        return false;
    }

	m_initizlied=true;
	return true;
}



}//end of namespace


