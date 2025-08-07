//#pragma once
#ifndef CUSG_DEF_PARSER_H
#define CUSG_DEF_PARSER_H

#include "def_struct.h"

namespace cusg
{

//-----------------------------------------------------------------------------
//This class parses definition files
//and also sets up and initialize the simulation

class DefinitionParser : public BaseParser
{
public:

    //constructor
    DefinitionParser();
    virtual ~DefinitionParser();

    //parsing the definition file
    bool parse(string filename);

	//
	virtual bool parse(istream& in, string folder="");
	
	//does nothing...
	virtual bool parse(block_raw_data& block){ return true; }
	
	//
	// access functions
	//
	block_raw_data * getBlockData(const string& name) { return  m_blocks[name];}
	map<string, block_raw_data *>& getBlockData() { return m_blocks; }

protected:
    
    typedef map<string, block_raw_data*(*)(block_raw_data& block)> BLOCK_TYPE_MAP;
    BLOCK_TYPE_MAP m_block_type_map;
    map<string, block_raw_data *> m_blocks;
};



template<typename T> block_raw_data* createBlockType(block_raw_data& block) { return new T(block); }

} //namespace cusg

#endif //CUSG_DEF_PARSER_H


