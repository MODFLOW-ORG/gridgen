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
    ~DefinitionParser() override;

    //parsing the definition file
    bool parse(string filename);

    //
    virtual bool parse(istream & in, string folder = "");

    //does nothing...
    virtual bool parse(DefinitionBlock & block) { return true; }

    //
    // access functions
    //
    DefinitionBlock * getBlockData(const string & name) { return m_blocks[name]; }
    map<string, DefinitionBlock *> & getBlockData() { return m_blocks; }

    //warnings for records that a block does not recognize
    const list<string> & getWarnings() const { return m_warnings; }

protected:
    using BLOCK_TYPE_MAP = map<string, DefinitionBlock * (*)(DefinitionBlock & block)>;
    BLOCK_TYPE_MAP m_block_type_map;
    map<string, DefinitionBlock *> m_blocks;
    list<string> m_warnings;
};

template <typename T> DefinitionBlock * createBlockType(DefinitionBlock & block)
{
    return new T(block);
}

} //namespace cusg

#endif //CUSG_DEF_PARSER_H
