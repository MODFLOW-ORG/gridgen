#pragma once
#ifndef CUSG_READER_TREESTRUCT_H
#define CUSG_READER_TREESTRUCT_H

#if defined(_WIN32) || defined(_WIN64)
#pragma warning(disable : 4244 4800 4267)
#endif

#include <string>
#include <vector>

using namespace std;

struct TreeStructNode
{
    TreeStructNode()
    {
        number = 0;
        layer = 0;
        row = 0;
        col = 0;
        code = "";
    }
    int number;
    int layer;
    int row;
    int col;
    string code;
};

//parse the parameters in the string node_str
void parse_node_str(string & node_str, TreeStructNode & node);

//this function is to parse the tree structure file into a quadtree instance
void parse_tree_struct_file(string file_name, vector<TreeStructNode> & tree,
                            string folder_path = "");

#endif