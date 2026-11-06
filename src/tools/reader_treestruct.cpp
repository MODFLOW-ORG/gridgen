#include "reader_treestruct.h"
#include "GridgenError.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
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
void parseNodeStr(string & nodeStr, TreeStructNode & node)
{
    trim(nodeStr);

    int number_end = nodeStr.find_first_of(',');
    string number_text = nodeStr.substr(0, number_end);

    int open_paren = nodeStr.find_first_of('(');
    int layer_end = nodeStr.find_first_of(',', open_paren);
    string layer_text = nodeStr.substr(open_paren + 1, layer_end - open_paren - 1);

    int row_end = nodeStr.find_first_of(',', layer_end + 1);
    string row_text = nodeStr.substr(layer_end + 1, row_end - layer_end - 1);
    int close_paren = nodeStr.find_first_of(')', row_end + 1);
    string col_text = nodeStr.substr(row_end + 1, close_paren - row_end - 1);

    if ((unsigned)(close_paren + 1) < nodeStr.size())
    {
        string location_text = nodeStr.substr(close_paren + 1, nodeStr.size() - close_paren - 1);
        for (char t : location_text)
        {
            if (t >= '1' && t <= '4') node.code.push_back(t);
        }
    }

    node.number = atoi(number_text.c_str());
    node.layer = atoi(layer_text.c_str());
    node.row = atoi(row_text.c_str());
    node.col = atoi(col_text.c_str());
}

void parseTreeStructFile(string fileName, vector<TreeStructNode> & tree, string folderPath)
{
    string filePath = folderPath;
    filePath += fileName;
    //replace the '\' with '/'
    int pathsize = filePath.size();
    for (int i = 0; i < pathsize; i++)
    {
        if (filePath[i] == '\\')
        {
            filePath[i] = '/';
        }
    }

    ifstream fin(filePath.c_str());
    if (!fin.good())
    {
        throw cusg::GridgenError("Cannot open tree structure file:" + filePath);
    } //end good

    /************************************************************************/
    //read the comment
    if (fin.eof()) return;

    string line;

    //Note: (3/25/2014) number of active nodes and number of connections are now stored in .nod file

    //read in total node number.
    int nodeNumber = -INT_MAX; //, active_nodeNumber=-INT_MAX, connectionNumber=-INT_MAX;

    while (nodeNumber == -INT_MAX) // || active_nodeNumber==-INT_MAX || connectionNumber==-INT_MAX)
    {
        getline(fin, line);
        while ((line.empty() || line[0] == '#') && !fin.eof())
        {
            getline(fin, line);
        }

        if (nodeNumber == -INT_MAX) nodeNumber = atoi(line.c_str());
    }

    for (int i = 0; i < nodeNumber; i++)
    {
        string lineStr;
        getline(fin, lineStr);
        if (!(line.empty() || line[0] == '#'))
        {
            tree.push_back(TreeStructNode());
            parseNodeStr(lineStr, tree.back());
        }
    } //end for i
}
