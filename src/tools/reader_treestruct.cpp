#include "reader_treestruct.h"
#include "GridgenError.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <limits.h>
using namespace std;

string & trim(string & s)
{
    if (s.empty())
    {
        return s;
    }

    string::iterator c;
    // Erase whitespace before the string
    for (c = s.begin(); c != s.end() && iswspace(*c++););
    s.erase(s.begin(), --c);

    // Erase whitespace after the string
    for (c = s.end(); c != s.begin() && iswspace(*--c););
    s.erase(++c, s.end());

    return s;
}

//parse the tree structure node
void parse_node_str(string & node_str, TreeStructNode & node)
{
    trim(node_str);

    int number_end = node_str.find_first_of(',');
    string number_text = node_str.substr(0, number_end);

    int open_paren = node_str.find_first_of('(');
    int layer_end = node_str.find_first_of(',', open_paren);
    string layer_text = node_str.substr(open_paren + 1, layer_end - open_paren - 1);

    int row_end = node_str.find_first_of(',', layer_end + 1);
    string row_text = node_str.substr(layer_end + 1, row_end - layer_end - 1);
    int close_paren = node_str.find_first_of(')', row_end + 1);
    string col_text = node_str.substr(row_end + 1, close_paren - row_end - 1);

    if ((unsigned)(close_paren + 1) < node_str.size())
    {
        string location_text = node_str.substr(close_paren + 1, node_str.size() - close_paren - 1);
        //any character but a space is kept, so a code that is not 1 to 4 is reported
        for (char t : location_text)
        {
            if (!isspace((unsigned char)t)) node.code.push_back(t);
        }
    }

    node.number = atoi(number_text.c_str());
    node.layer = atoi(layer_text.c_str());
    node.row = atoi(row_text.c_str());
    node.col = atoi(col_text.c_str());
}

void parse_tree_struct_file(string file_name, vector<TreeStructNode> & tree, string folder_path)
{
    string file_path = folder_path;
    file_path += file_name;
    //replace the '\' with '/'
    int pathsize = file_path.size();
    for (int i = 0; i < pathsize; i++)
    {
        if (file_path[i] == '\\')
        {
            file_path[i] = '/';
        }
    }

    ifstream fin(file_path.c_str());
    if (!fin.good())
    {
        throw cusg::GridgenError("cannot open structure file " + file_path);
    } //end good

    /************************************************************************/
    //read the comment
    if (fin.eof()) return;

    string line;

    //Note: (3/25/2014) number of active nodes and number of connections are now stored in .nod file

    //read in total node number.
    int node_number = -INT_MAX; //, active_nodeNumber=-INT_MAX, connectionNumber=-INT_MAX;

    while (node_number == -INT_MAX) // || active_nodeNumber==-INT_MAX || connectionNumber==-INT_MAX)
    {
        getline(fin, line);
        while ((line.empty() || line[0] == '#') && !fin.eof())
        {
            getline(fin, line);
        }

        if (node_number == -INT_MAX) node_number = atoi(line.c_str());
    }

    for (int i = 0; i < node_number; i++)
    {
        string line_str;
        getline(fin, line_str);
        if (!(line.empty() || line[0] == '#'))
        {
            tree.push_back(TreeStructNode());
            parse_node_str(line_str, tree.back());
        }
    } //end for i
}
