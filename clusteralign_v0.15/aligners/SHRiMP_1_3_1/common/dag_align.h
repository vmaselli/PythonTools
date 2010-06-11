/*	$Id: dag_align.h 245 2008-06-06 18:24:28Z rumble $	*/

#ifndef ALIGN_H_
#define ALIGN_H_

#include <string> 

typedef enum {MATCH, LETTERGAP, GAPGAP, MISMATCH} pair_type;
//typedef enum {FORWARD, REVERSE} Direction;

//union Letter;
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <climits>
using namespace std; 


class Edge;
class Column;
class Node;
class AlignNode;
class AlignEdge;
class CrossProduct;

class kmersInfo
{
public:
	vector <set <string> > kmers;
	kmersInfo(int maxsize);
	void mergeinto (kmersInfo kmers);
};

class Alignment
{
public:
	int seqstart; //staring position in the genome window. Starts with zero 
	int seqend;
	string seq; //the genome sequence as appears in the alignment
	string read1;
	string read2;
	Alignment();
};

class Node
{
public:
	vector <Edge*> pred_list; //neighbor
	vector <Edge*> succ_list; //neighbor 
	const int index; //index in the node_list vector 
	kmersInfo *aux;	
	Node(int index);
	void AddSelfLoop(int column_length, char );
	
};

ostream& operator << (ostream& os, Node& n);

class Graph
{
public:
	const int column_length;
	int num_nodes;
	vector <Node*> node_list;
	int maxkmersize;
	Graph (string s);
	Graph (CrossProduct & g, int epsilon);
	void AddSelfLoops();
	void mergeinto (Node *source, Node *dest, char letter, kmersInfo & kmers);
	kmersInfo *getkmers(int size);	
	void printkmers(int size);	
	~Graph();
};


ostream& operator << (ostream& os, Graph& graph);

class AlignNode
{
public:	
	vector <AlignEdge*> pred_list; //neighbor
	vector <AlignEdge*> succ_list; //neighbor 
	int x,y; //index in the CrossProduct array?
	int fscore, bscore;
	AlignNode* dijkstra_parent;
	AlignEdge* dijkstra_parent_edge;
	Node* gnode; /* This is used to point to the location of the corresponding 
	Graph node in the constructor Graph (CrossProduct & g, int epsilon).
	Much better would be to use a map from (i,j) to Node*. Should be redone */ 
	
	AlignNode();
	AlignNode (int x, int y);
};

class SmallAlignNode
{
public:	
	int fscore;
	SmallAlignNode* parent;
	SmallAlignNode();
};

ostream& operator << (ostream& os, AlignNode& n);

class SmallCrossProduct
{
public:
	Graph & a;
	Graph & b;
	vector <vector <SmallAlignNode> >   matrix;
	
	/* Used for a crude iterator over pairs of nodes. 
    The iterators are initialized in the constructor. They are initialized
    every time iteration through the nodes completes. Can be used 
	for a single iteration of pairs in a time */
	vector<Node*>::iterator ai;
	vector<Node*>::iterator bi;
	inline bool VertexPairIterator(Node ** nap, Node ** nbp, SmallAlignNode ** san);
	
	/*Used for iteration over pairs of predecessor edges */
	vector<Edge*>::iterator predai, predbi;
	inline bool EdgePredPairIterator(Node * nap, Node * nbp, Edge ** eapp, Edge ** ebpp, SmallAlignNode **sanpp);
	
	const int column_length;
	
	SmallCrossProduct (Graph & a, Graph & b);
	int DijkstraForward (SmallAlignNode ** pbestend, Node ** bestendapp, Node ** bestendbpp);
	void smallinitscore();
	int PrintBestPath(SmallAlignNode * bestend, Node * bestendap, Node * bestendbp);
	bool EdgesBetween(SmallAlignNode * sourcep, SmallAlignNode * destp, 
				Node * destap, Node * destbp, Edge ** eapp, Edge ** ebpp);
	Alignment * GetBestPath(SmallAlignNode * bestend, Node * bestendap, Node * bestendbp);
  //  ~SmallCrossProduct();
};

class CrossProduct
{
public:
	Graph & a;
	Graph & b;
	vector <vector <AlignNode*> >   matrix;
	const int column_length;		
	CrossProduct (Graph & a, Graph & b);
	int DijkstraForward (bool global,  AlignNode ** pbestend);
	int DijkstraBackward (bool global);//,  AlignNode * bestend);
	void initscore(bool global);
	int PrintBestPath (AlignNode* bestend);
    ~CrossProduct();
};

ostream& operator << (ostream& os, CrossProduct& g);
ostream& operator << (ostream& os, Edge& edge);

class Column
{
public:
	string s;
		
	
	
	static  int CC;       //Score of aligning two equal letters
	static  int C_ ;
	static  int CnotC; /* Score of aligning two not equal letters - 
	                             we disallow alignment of non matching characters */

	/* DAG to Genome Alignment */
	static int S_MATCH;
	static int S_SNP;
	static int S_HALF_MATCH;
	static int S_NEITHER_MATCH;
	static int S_MATCH_DELETION;
	static int S_MISMATCH_DELETION;
	static int S_ERROR_INSERTION;
	
	/* Score of the alignment of two letters.*/
	 int get2score (char c1, char c2);
	

	/* Finds the most frequent character in a column. Return arbitrary 
	 * 	if they are equally likely */
	const int score;
	char consensus (int Acount, int Ccount, int Gcount,
			int Tcount, int _count, Column *col2p)	;
	char consensus();
	Column (string s);
	Column (Column &c1, Column &c2);
	Column (const Column &c1);
	int getscore(Column* col2=NULL);
	static void setscore(int read_match, int read_gap, int read_mismatch, int dag_match,
	    int dag_snp, int dag_half_match, int dag_neither_match, int dag_match_deletion,
	    int dag_mismatch_deletion, int dag_error_insertion);
       	int get3score(Column* col2);
	static pair_type get_pair_type(Column * col2, char * letter); 	
};

ostream& operator << (ostream& os, Column& c);

class Edge
{
	public:
	Node & source;
	Node & dest;
	Column c;
	Edge (Node & source, Node & dest, string s);
	Edge (Node & src, Node & dest, Column & c);
	static int score(Edge * eap, Edge * ebp);
};




class AlignEdge
{
public:
	AlignNode & src, & dest;
	Edge* e1, *e2;
	Column col;
	int score;
	AlignEdge (AlignNode & src, AlignNode & dest,  Edge * e1, Edge * e2);
	AlignEdge( const AlignEdge& lhs);
};

ostream& operator << (ostream& os, AlignEdge& e);
ostream& operator << (ostream& os, kmersInfo& e);

#endif /* ALIGN_H_ */
