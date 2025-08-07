#include "reader_treestruct.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <limits.h>
using namespace std;

string& trim(string& s) {
	  if (s.empty()) {
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
void parseNodeStr(string& nodeStr, TreeStructNode& node)
{
	trim(nodeStr);

	int comIdx1 = nodeStr.find_first_of(',');
	string numStr = nodeStr.substr(0,comIdx1);

	int comIdx2 = nodeStr.find_first_of('(');
	int comIdx3 = nodeStr.find_first_of(',',comIdx2);
	string layerstr = nodeStr.substr(comIdx2 + 1, comIdx3 - comIdx2 - 1);

	int comIdx4 = nodeStr.find_first_of(',', comIdx3 + 1);
	string rowstr = nodeStr.substr(comIdx3 + 1, comIdx4 - comIdx3 - 1);
	int comIdx5 = nodeStr.find_first_of(')', comIdx4 + 1);
	string colstr = nodeStr.substr(comIdx4 + 1,comIdx5 - comIdx4  - 1);

	if((unsigned)(comIdx5 + 1) < nodeStr.size())
	{
		string tmpCandStr = nodeStr.substr(comIdx5 + 1,nodeStr.size() - comIdx5 - 1);
		for(string::iterator site = tmpCandStr.begin(); site != tmpCandStr.end(); ++site)
		{
			char t = *site;
			if(t>='1'&&t<='4')
				node.code.push_back(t);
		}
	}

	node.number = atoi(numStr.c_str());
	node.layer = atoi(layerstr.c_str());
	node.row = atoi(rowstr.c_str());
	node.col = atoi(colstr.c_str());
}

void parseTreeStructFile(string fileName, vector<TreeStructNode> & tree,string folderPath)
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
	if(!fin.good()){
		cerr<<"! Error: Cannot open tree structure file:" <<filePath<<endl;
		exit(1);
		return;
		//return false;
	}//end good

	/************************************************************************/
	//read the comment
	if(fin.eof())
		return;

	string line;
	
	//Note: (3/25/2014) number of active nodes and number of connections are now stored in .nod file

	//read in total node number.
	int nodeNumber=-INT_MAX;//, active_nodeNumber=-INT_MAX, connectionNumber=-INT_MAX;

	while(nodeNumber==-INT_MAX)// || active_nodeNumber==-INT_MAX || connectionNumber==-INT_MAX)
	{
		getline(fin, line);
		while( (line.empty() || line[0] == '#') && !fin.eof())
		{
			getline(fin, line);
		}

		if(nodeNumber==-INT_MAX) nodeNumber = atoi(line.c_str());
		//else if(active_nodeNumber==-INT_MAX) active_nodeNumber = atoi(line.c_str());
		//else if(connectionNumber==-INT_MAX) connectionNumber = atoi(line.c_str());
	}

	for (int i = 0;i < nodeNumber; i++)
	{
		string lineStr;
		getline(fin, lineStr);
		if(!(line.empty() || line[0] == '#' ) )
		{
			tree.push_back(TreeStructNode());
			parseNodeStr(lineStr,tree.back());
		}
	}//end for i
}
