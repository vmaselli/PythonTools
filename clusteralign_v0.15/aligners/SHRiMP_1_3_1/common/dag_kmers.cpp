/*	$Id: dag_kmers.cpp 235 2008-05-15 16:32:08Z rumble $	*/

#include <string> 
#include "dag_align.h"
//typedef enum {FORWARD, REVERSE} Direction;

//union Letter;
#include <iostream>
#include <vector>
#include <string>
#include <climits>
#include <fstream>

using namespace std;


kmersInfo::kmersInfo(int maxsize)
{
	kmers.resize(maxsize,  set <string>() );
}



/* Assumes dest is allocated! */
void Graph::mergeinto (Node *source, Node *dest,  char character, kmersInfo & kmers)
{
	int count = 0;
	string letter (1, character);
	if (dest->aux == NULL)
	{
		 dest->aux = new kmersInfo (maxkmersize);
	}
	set <string>::iterator setit;
	((dest->aux->kmers).at(0)).insert(letter) ;
	kmers.kmers.at(count).insert (letter);
	for(count = 0;count <= maxkmersize -1; count++ )
	{
		//fromit = &((source->aux->kmers).at(count));
		//toit=&((dest->aux->kmers).at(count));
		for(setit= (source->aux->kmers).at(count).begin(); setit != (source->aux->kmers).at(count).end(); setit++)
		{
			kmers.kmers.at(count).insert (*setit);
			if(count<=maxkmersize -2)
				((dest->aux->kmers).at(count+1)).insert(*setit+letter) ;
		}
		
	}
}

kmersInfo *Graph::getkmers(int size)
{
	maxkmersize = size;
	vector<Node*>::iterator ni;
	vector<Edge*>::iterator ei;
	ni = node_list.begin();
	
	kmersInfo * kmers = new kmersInfo (maxkmersize);
	node_list[0]->aux = new kmersInfo (maxkmersize);
	
	for (; ni!= node_list.end(); ni++)
	{
		for ( ei = (**ni).succ_list.begin(); ei != (**ni).succ_list.end(); ei++)
		{
			char letter = ((**ei).c).consensus();
			/*if (&((**ei).dest)==&((**ei).source))
				cout << "voila\n"; */
			mergeinto((*ni), &((**ei).dest), letter, *kmers);
		}
		delete (**ni).aux;
	}
	
	return (kmers);
}

void Graph::printkmers(int size)
{
	kmersInfo *kmers = getkmers(size); 
	cout << *kmers;
	delete kmers;
}

ostream& operator << (ostream& os, kmersInfo& info)
{
	set <string>::iterator setit;
	for(unsigned count = 0; count <= info.kmers.size()-1; count++ )
		{
			//fromit = &((source->aux->kmers).at(count));
			//toit=&((dest->aux->kmers).at(count));
			for(setit= (info.kmers).at(count).begin(); setit != (info.kmers).at(count).end(); setit++)
			{
				os << (*setit)<<endl;
			}
		}
	return os;
}
