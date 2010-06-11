/*	$Id: dag_align.cpp 292 2008-11-18 20:50:53Z rumble $	*/

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


int MINSCORE = -1000;

/* Read to Read Alignment */
int Column::CC = -1000;		//2;
int Column::C_ = -1000;		//-1;
int Column::CnotC = -1000;	//-100;

/* DAG to Genome Alignment */
int Column::S_MATCH 		= -1000;	//11;
int Column::S_SNP	  	= -1000;	//-10;
int Column::S_HALF_MATCH	= -1000;	//4;
int Column::S_NEITHER_MATCH	= -1000;	//-5;
int Column::S_MATCH_DELETION	= -1000;	//5;
int Column::S_MISMATCH_DELETION	= -1000;	//-6;
int Column::S_ERROR_INSERTION	= -1000;	//-6;

void Column::setscore(int read_match, int read_gap, int read_mismatch, int dag_match,
    int dag_snp, int dag_half_match, int dag_neither_match, int dag_match_deletion,
    int dag_mismatch_deletion, int dag_error_insertion)
{
	Column::CC = read_match;
	Column::C_ = read_gap;
	Column::CnotC = read_mismatch;

	Column::S_MATCH			= dag_match;
	Column::S_SNP			= dag_snp;
	Column::S_HALF_MATCH		= dag_half_match;
	Column::S_NEITHER_MATCH		= dag_neither_match;
	Column::S_MATCH_DELETION	= dag_match_deletion;
	Column::S_MISMATCH_DELETION	= dag_mismatch_deletion;
	Column::S_ERROR_INSERTION	= dag_error_insertion;
}

int Column::get2score (char c1, char c2)
	{
		if (c1 == '-' || c2=='-')
			return  0; // no penalty for gaps at ends
		if (c1 == '_' || c2=='_')
			return Column::C_;
		if (c1 == c2) 
			return Column::CC;
		        
		return Column::CnotC; //CnotC is very high to disallow mismatches   
	}


/*
char Column::consensus ()
	{
		int Acount=0, Ccount=0, Gcount=0, Tcount=0, _count=0;
		
		string::const_iterator siter = s.begin();
		
		while ( siter != s.end() ) 
		{
				   switch (*siter)
				   {
				   case 'C':
					   Ccount++; break;
				   case 'T':
					   Tcount++; break;
				   case 'G':
					   Gcount++; break;
				   case 'A':
					   Acount++; break;
				   case '_':
                                   case '-':
					   _count++; break;
				   default: cout << "Seeing very strange character"<<*siter<<". Help!"<< endl;	   
				     
				   }
				   siter++;
		}
		int max = Acount;  		char result = 'A';
		if (Ccount > max) {max = Ccount;  result = 'C';}
		if (Gcount > max) {max = Gcount;  result = 'G';}
		if (Tcount > max) {max = Tcount;  result = 'T';}
		if (_count > max) {max = _count;  result = '_';}
        return result;
	} */
char Column::consensus()
{
	return Column::consensus(0,0,0,0,0,NULL);
}


char Column::consensus (int Acount=0, int Ccount=0, int Gcount=0, int Tcount=0, int _count=0, Column* col2 = NULL)
	{
		string::const_iterator siter = s.begin();
		
		while ( siter != s.end() ) 
		{
				   switch (*siter)
				   {
				   case 'C':
					   Ccount++; break;
				   case 'T':
					   Tcount++; break;
				   case 'G':
					   Gcount++; break;
				   case 'A':
					   Acount++; break;
				   case '_':
                                   case '-':
					   _count++; break;
				   default: cout << "Seeing very strange character"<<*siter<<". Help!"<< endl;	   
				     
				   }
				   siter++;
		}
		int max = Acount;  		char result = 'A';
		if (Ccount > max) {max = Ccount;  result = 'C';}
		if (Gcount > max) {max = Gcount;  result = 'G';}
		if (Tcount > max) {max = Tcount;  result = 'T';}
		if (_count > max) {max = _count;  result = '_';}
        
		if (col2 == NULL)
			return result;
		else return col2->consensus (Acount, Ccount, Gcount, Tcount, _count, NULL );
	}


Column::Column (string s):s(s),score(getscore())
{

}
Column::Column (Column &c1, Column &c2):s(c1.s + c2.s),score(getscore())
{

}

Column::Column (const Column &c1):s(c1.s), score(c1.score)
{

}

int Column::getscore(Column* col2)
{
	int result = 0;
	if (col2 == NULL && s.length()==2)
			return get2score(s[0],s[1]);
	string::const_iterator siter = s.begin();
	//Rewrite so that reallocation of a column is not needed to compute the score
	char cons = consensus (0,0,0,0,0,col2);
	while ( siter != s.end() ) 
	{
	   result += get2score (cons, *siter); //get2score should be defined!!!!!!!!!!!! 
	   siter++;        
	}
	
	if (col2 != NULL)
	{
		siter = (col2->s).begin();
		while ( siter != (col2->s).end() ) 
			{
			   result += get2score (cons, *siter); //get2score should be defined!!!!!!!!!!!! 
			   siter++;        
			}
	}
	return result;
}

pair_type Column::get_pair_type(Column * col2, char * letter)
{
	
	if ((col2->s).length() != 2)
	{
		cout <<"get_pair_type: bad invocation\n";
		return MISMATCH;
	}
	char l0 = (col2->s)[0], l1 = (col2->s)[1];
	if ((l0 == '_' ||l0 == '-') && (l1 == '_' ||l1 == '-'))
	{
		*letter = l0;
		return GAPGAP;
	}
	
	if (l0 == '_' ||l0 == '-') 
	{
		*letter = l1;
		return LETTERGAP;
	}

	if (l1 == '_' ||l1 == '-')
	{
		*letter = l0;
		return LETTERGAP;
	}
	if (l0 == l1)
	{
		*letter = l0;
		return MATCH;
	}
	
	if (l0 != l1)
	{
		*letter = l0;
		return MISMATCH;
	}

	/* shut up the compiler */
	cout << "IMPOSSIBLE!! FELL THROUGH!!" << endl;
	return MISMATCH;
}
 
int Column::get3score(Column* col2)
{
	pair_type p;
	
	/* not generic here */
	if ((col2->s).length() != 2 || s.length()!=1 )
	{
		cout << "get3score: Expecting to match a genome vs the read pair\n";
		return MINSCORE;
	}
	
	char letter; //the letter representing the read pair
	
	char gen = s[0]; //the letter in the genome. Can be A, C, G, T, _, - (- becaus addselfloop does not care who you are)
	if (gen == '-')
	{
		gen = '_'; // adjusting for the "special" gaps
	}
	
	p = Column::get_pair_type(col2, &letter);
	
	
	switch (p)
	{
	case MATCH:
		if (gen == letter) // A AA
			return S_MATCH;
		else
			return S_SNP; //B AA
		
	case LETTERGAP:
		if (gen == '_')
		{
			return S_ERROR_INSERTION; // _ _A
		}
		
		if (gen == letter) //A _A
		{
			return S_MATCH_DELETION;
		}
		else // A _C
		{
			return S_MISMATCH_DELETION;
		}
		
	case GAPGAP:
		return S_NEITHER_MATCH; //A _ _    or    _  _ _
	case MISMATCH:
		if ((gen == (col2->s)[0])||(gen == (col2->s)[1]))
			return S_HALF_MATCH;                //A 
		else
			return S_NEITHER_MATCH;
		
	default:
		return MINSCORE;
	}

}



ostream& operator << (ostream& os, Column& c)
{
	os << c.s;
	return os;
}

   
Edge::Edge (Node & source, Node & dest, string s):source(source), dest(dest), c(s)
{
	(source.succ_list).push_back(this);	
	(dest.pred_list).push_back(this);	
}



Edge::Edge (Node & source, Node & dest, Column & c):source(source), dest(dest), c(c)
{
	(source.succ_list).push_back(this);	
	(dest.pred_list).push_back(this);	
}

ostream& operator << (ostream& os, Edge& edge)
{
	os << "Edge between " << edge.source << " and " << edge.dest  <<" is " << edge.c;
	return os;
}

ostream& operator << (ostream& os, Node& n)
{
	os << n.index;
	return os;
}

Graph::Graph (string s):column_length(1)
{
	Node* prev = new Node(0);
	num_nodes=1;
	string::iterator i = s.begin();
	node_list.push_back (prev);
	while (i != s.end())
	{
		Node* next = new Node(num_nodes); //num_nodes initializes index field
		node_list.push_back (next);
		new Edge(*prev, *next, string(1,*i));
		prev = next;
		i++;
		num_nodes++;
	}
} 

Graph::~Graph()
{
	vector<Node*>::iterator ni;
	vector<Edge*>::iterator ei;
	for (ni = node_list.begin(); ni!= node_list.end(); ni++)
	{
		for ( ei = (**ni).succ_list.begin(); ei != (**ni).succ_list.end(); ei++)
		{
			delete *ei;
		}
		
	}
	
	for (ni = node_list.begin(); ni!= node_list.end(); ni++)
	{
		delete *ni;
	}
}


Graph::Graph(CrossProduct & g, int epsilon):column_length(g.column_length)
{
//	cout << "Hello from Graph::Graph(CrossProduct & g, int epsilon)" << endl;
	num_nodes=0;
	vector<Node*>::iterator ai, bi;
	vector<AlignEdge*>::const_iterator ei;
	
	AlignNode * align_dest;
	
	/* g[0][0] should not have any incoming edges so it would not be added twice
	 * but it must be in the graph built anyway*/
	
	Node* dest = new Node(num_nodes); 
	node_list.push_back (dest);
	Node* src;

	bool global = true;
	int bestforward, bestbackward, unneeded_var_checked_nodes=0;
	AlignNode * bestend; //Does not matter for global alignment
	g.initscore( global);
	bestforward = g.DijkstraForward(global, &bestend);
	bestbackward = g.DijkstraBackward(global);//, bestend);
	if (bestforward != bestbackward) {
		cout << "Help! Something very wrong:"<<"forward="<<bestforward<<" backward="<<bestbackward<<endl;
	}
	else ;//cout << "bestforward = " << bestforward << endl; 

	g.matrix[0][0]->gnode = dest;
	for (ai=g.a.node_list.begin(); ai!= g.a.node_list.end(); ai++) {
		for (bi=g.b.node_list.begin(); bi!= g.b.node_list.end(); bi++) {
			align_dest = g.matrix[(**ai).index][(**bi).index];
			ei = (*align_dest).pred_list.begin();
			unneeded_var_checked_nodes++;
			if (ei == (*align_dest).pred_list.end())
				continue;
			
			
			//cout << "\n\nWorking on g.matrix["<< (**ai).index<<"]["<<(**bi).index<<"] here which is "<<*align_dest<<endl;
			bool dest_created = false;
			/* dest_created is false for matrix[0][0] though the node was created. It
			 * does not matter since  matrix[0][0] has no predecessors */
			for (; ei != (*align_dest).pred_list.end(); ei++) {
				//cout << "working on edge "<< **ei<<" connecting to a predecessor "<<*align_dest << endl;
				if (((**ei).src).fscore + (**ei).score + ((**ei).dest).bscore >= bestforward - epsilon &&
				     &((**ei).src) != &((**ei).dest))
				{
				//	cout << "Seeing good edge with score: "<<
						//((**ei).src).fscore + (**ei).score + ((**ei).dest).bscore <<endl;
					//cout << "src's fscore = "<< ((**ei).src).fscore<<" dest bscore="<<((**ei).dest).bscore;
					//cout <<" edge score = "<<  (**ei).score<<endl;
					/* Good edge. If the destination not created yet - create */
					if (!dest_created) {
						num_nodes++;
						dest_created = true;
						//cout << "Creating node with index "<<num_nodes <<endl;
						dest = new Node(num_nodes);
						
						node_list.push_back (dest);
						g.matrix[(**ai).index][(**bi).index]->gnode = dest;
					}
					src = ((**ei).src).gnode;
					//cout << "src index="<<src->index<<";dest index="<<dest->index<<endl;
					// cout << "src.x="<<((**ei).src).x<<";dest.x="<< ((**ei).dest).x<<endl;
					new Edge(*src, *dest, (*ei)->col);
					//cout << "Edge between "<< *src << " and "<< *dest <<" was created\n"; 
				}
			}

		}
	}
	num_nodes++; //It was initialized to zero!!!!!
}		


ostream& operator << (ostream& os, AlignEdge& edge)
{
	os << "edge between " << edge.src << " and " << edge.dest  <<" is " << edge.col<<" of score "<<edge.score<<endl;
		return os;
}

ostream& operator << (ostream& os, Graph& graph) 
{
	vector<Node*>::iterator ni;
	vector <Edge*>::iterator ei;
	for(ni = graph.node_list.begin(); ni<graph.node_list.end();ni++)
	{
		cout <<"Printing successors of "<< **ni<<endl;
		for ( ei = (**ni).succ_list.begin(); ei != (**ni).succ_list.end(); ei++)
		{
			cout << (**ei)<<endl;
		} 
		cout <<"end of successors of "<< **ni<<endl;
		
	}
	return os;
}


CrossProduct::~CrossProduct() {

	AlignNode *dest;
	vector<AlignEdge*>::const_iterator ei;
	//initscore(forward, global);  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	vector<Node*>::iterator ai;
	vector<Node*>::iterator bi;

	for (ai=a.node_list.begin(); ai!= a.node_list.end(); ai++) {
		for (bi=b.node_list.begin(); bi!= b.node_list.end(); bi++) {
			dest = matrix[(**ai).index][(**bi).index];
			
		
			ei = (*dest).pred_list.begin();
			if (ei!= (*dest).pred_list.end())
			 {
				for (; ei != (*dest).pred_list.end(); ei++) {
					delete *ei;
			}
			}
		}
	}

	for (ai=a.node_list.begin(); ai!= a.node_list.end(); ai++) {
			for (bi=b.node_list.begin(); bi!= b.node_list.end(); bi++) {
				delete matrix[(**ai).index][(**bi).index];
			}
	}
}

ostream& operator <<(ostream& os, CrossProduct& g) {

	AlignNode *dest;
	vector<AlignEdge*>::const_iterator ei;
	//initscore(forward, global);  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	vector<Node*>::iterator ai;
	vector<Node*>::iterator bi;

	for (ai=g.a.node_list.begin(); ai!= g.a.node_list.end(); ai++) {
		for (bi=g.b.node_list.begin(); bi!= g.b.node_list.end(); bi++) {
			dest = g.matrix[(**ai).index][(**bi).index];
			os <<endl<<"Preds of "<<*dest ;
			os <<" fscore="<<(*dest).fscore <<" bscore="<<(*dest).bscore<< " parent= ";
			if ((*dest).dijkstra_parent) {
				os<< "("<<(*dest).dijkstra_parent->x<<";"<<(*dest).dijkstra_parent->y<<")" << endl;
			} else {
				os<<"NULL" <<endl;//<<";"<<n.dijkstra_parent->y<<")";*/; 
			}
 
			ei = (*dest).pred_list.begin();
			if (ei == (*dest).pred_list.end())
				os<< "No predecessors";
			else {
				for (; ei != (*dest).pred_list.end(); ei++) {
					os << **ei;
				}
			}
		}
	}

	return os;

}


Node::Node(int index):index(index), aux(NULL)
{
}

void Node::AddSelfLoop(int column_length, char c)
{
	new Edge(*this, *this, string(column_length,c));
}

void Graph::AddSelfLoops()
		{
			vector<Node*>::iterator ai;
			
			for (ai=node_list.begin(); *ai!= node_list.back(); ai++) {
				(**ai).AddSelfLoop(column_length,'_');
			}
			
			(**ai).AddSelfLoop(column_length,'-');
		}


SmallAlignNode::SmallAlignNode():fscore(0), parent(NULL){}
AlignNode::AlignNode():x(0),y(0),fscore(0),bscore(0){}
AlignNode::AlignNode (int x, int y):x(x), y(y), fscore(0),bscore(0){}
ostream& operator << (ostream& os, AlignNode& n)
{
	os <<"("<<n.x<<";"<<n.y<<")"  ;
	return os;
}

AlignEdge::AlignEdge (AlignNode & src, AlignNode & dest,  Edge * e1, Edge * e2):src(src),dest(dest), e1(e1), e2(e2), col(e1->c,e2->c)
	{
		score = col.score ;
	}


AlignEdge::AlignEdge( const AlignEdge& lhs) : src(lhs.src),dest(lhs.dest), e1(lhs.e1), e2(lhs.e2), col(lhs.e1->c,lhs.e2->c) 
{	
}

/*
AlignEdge::AlignEdge & operator=( const AlignEdge& lhs ) {
		if( &lhs == this ) return *this;
		
		src = lhs.src;
		dest = lhs.dest;
		e1 = lhs.e1;
		e2 = lhs.e2;
		
	}
*/



CrossProduct::CrossProduct (Graph & a, Graph & b): a(a),b(b),column_length(a.column_length + b.column_length) 
{
	int asize = a.num_nodes;
	int bsize = b.num_nodes;
	matrix.resize(asize,  vector<AlignNode*>(bsize));
	vector<Node*>::iterator ai;
	vector<Node*>::iterator bi;
	vector<Edge*>::iterator sai, sbi;
	
	a.AddSelfLoops();
	b.AddSelfLoops();
//	cout << "Welcome from Cross constructor; a="<<asize<<" b="<<bsize<<endl;
	//cout << "a=\n"<<a<<"\nb=\n"<<b;
	for (ai=a.node_list.begin(); ai!= a.node_list.end(); ai++)
	{
		for (bi=b.node_list.begin(); bi!= b.node_list.end(); bi++)
		{
			AlignNode* dest = new AlignNode((**ai).index,(**bi).index);
			matrix[(**ai).index][(**bi).index] = dest;
			for( sai = ((**ai).pred_list).begin(); sai != ((**ai).pred_list).end(); ++sai) {  //scc1i - iterator over the vector of successors of node *g1i
				for (sbi = (**bi).pred_list.begin(); sbi != ((**bi).pred_list).end(); ++sbi)
				{
					AlignNode* src = matrix[((**sai).source).index][((**sbi).source).index] ; 
					AlignEdge* pe= new AlignEdge (*src, *dest, *sai, *sbi);
					(src->succ_list).push_back(pe);
					(dest->pred_list).push_back(pe);
				}
			}
		}
	}
}



int CrossProduct::DijkstraForward (bool global,  AlignNode ** pbestend)
{
	
	AlignNode *dest, *src;
	vector<AlignEdge*>::const_iterator ei;
	vector<Node*>::iterator ai;
	vector<Node*>::iterator bi;
	int maxscore = MINSCORE;
	AlignNode * bestend = NULL;

	/* shut up, icc */
	(void)global;
	
	for (ai=a.node_list.begin(); ai!= a.node_list.end(); ai++) 
	{
		for (bi=b.node_list.begin(); bi!= b.node_list.end(); bi++)
		{
			dest = matrix[(**ai).index][(**bi).index] ;
			for ( ei = (*dest).pred_list.begin(); ei != (*dest).pred_list.end(); ei++)
			{
				src = &((**ei).src);
				if (src->fscore > MINSCORE && (**ei).score+src->fscore>dest->fscore)
				{
					dest->fscore = (**ei).score+src->fscore;
					dest->dijkstra_parent = src; 
					dest->dijkstra_parent_edge = *ei;
					//cout << "dijkstra_parent="<<*src<<endl;
					if (dest->fscore > maxscore)
					{
						maxscore = dest->fscore;
						bestend = dest;
						//cout << "assigned maxscore="<<maxscore<<endl<<"bestend="<<*bestend<<endl;
						
					}
				}
			}
		}
	}
	*pbestend =bestend ;
//	cout << "Finally, maxscore="<<maxscore<<"bestend="<<*bestend<<endl;
	return matrix[a.node_list.back()->index][b.node_list.back()->index]->fscore;	
}

int SmallCrossProduct::DijkstraForward (SmallAlignNode ** pbestend, Node ** bestendapp, Node ** bestendbpp)
{
	
	SmallAlignNode *destp, *srcp;
	Edge * eap, * ebp;
	Node* nap;
	Node* nbp;
	
	int maxscore = MINSCORE;
	SmallAlignNode * bestend = NULL;
	while (VertexPairIterator(&nap,  &nbp,  &destp))
	{
		predai = (*nap).pred_list.begin();
		predbi = (*nbp).pred_list.begin();
		while(EdgePredPairIterator( nap, nbp, &eap, &ebp, &srcp))
		{
			int escore =  (eap->c).get3score(&(ebp->c));
			if (srcp->fscore > MINSCORE && escore + srcp->fscore  > destp->fscore)
			{
				destp->fscore = srcp->fscore + escore;
				destp->parent = srcp;
				if (destp->fscore > maxscore)
				{
					maxscore = destp->fscore;
					bestend = destp;
					*bestendapp = nap;
					*bestendbpp = nbp;
					//cout << maxscore<<endl;
				}
				
			}
		}
	}
	
	*pbestend = bestend ;
	return maxscore;
}



int CrossProduct::PrintBestPath(AlignNode * bestend)
{
	int len = 0;
	if (bestend == NULL)
		cout << "bestend == NULL!!!!!!!!!!!!!!!";
	else ; //cout << "Outputing the best alignment\n";
	
	for (AlignNode* current = bestend; current->dijkstra_parent_edge != NULL; current = current->dijkstra_parent)
	{
		cout<<(current->dijkstra_parent_edge)->col<<endl;
	}
	return len;
}

int SmallCrossProduct::PrintBestPath(SmallAlignNode * bestend, Node * bestendap, Node * bestendbp)
{
	int len = 0;
	Edge *ea = NULL, *eb = NULL;
	if (bestend == NULL)
		cout << "smallbestend == NULL!!!!!!!!!!!!!!!";
	else ; //cout << "Outputing the best alignment\n";
	
	for (SmallAlignNode* current = bestend; current->parent != NULL;
		current = current->parent, bestendap = &(ea->source), bestendbp = &(eb->source))
	{
		if (!EdgesBetween(current->parent, current, bestendap, bestendbp, &ea, &eb ))
			{cout << "EdgeBetween couldn't find edge\n";}
		else
			{
				cout<< (ea->c).s<<" "<<(eb->c).s<<endl;
				
			}
		//cout << len++;
	}
	cout << len;
	return len;
}

Alignment::Alignment()
{
	seqstart=0;
	seqend=0;
	string seq="";
	string read1="";
	string read2="";
}

Alignment* SmallCrossProduct::GetBestPath(SmallAlignNode * bestend, Node * bestendap, Node * bestendbp)
{

	Edge *ea = NULL, *eb = NULL;
	if (bestend == NULL)
		cout << "smallbestend == NULL!!!!!!!!!!!!!!!";
	else ; //cout << "Outputing the best alignment\n";
	
	Alignment *align = new Alignment(); 
	align->seqend = bestendap->index - 1;
	
	for (SmallAlignNode* current = bestend; current->parent != NULL;
		current = current->parent, bestendap = &(ea->source), bestendbp = &(eb->source))
	{
		if (!EdgesBetween(current->parent, current, bestendap, bestendbp, &ea, &eb ))
			{cout << "EdgeBetween couldn't find edge\n";}
		else
			{
				if ((ea->c).s.length()!=1)
				{
					cout << "The genome must be the first argument to Dijkstra, not reads' grapn\n";
				}
				else
				{
					align->seq = (ea->c).s+align->seq;
				}
				
				
				if ((eb->c).s.length()!=2)
				{
					cout << "Working with read pairs as of now\n";
				}
				else
				{
					align->read1 = string(1,((eb->c).s)[0])+align->read1;			
					align->read2 = string(1,((eb->c).s)[1])+align->read2;
				}
			}
		}
	align->seqstart = bestendap->index;
	return align;
}

bool SmallCrossProduct::EdgesBetween(SmallAlignNode * sourcep, SmallAlignNode * destp, 
		Node * destap, Node * destbp, Edge ** eapp, Edge ** ebpp)
{
	SmallAlignNode * temp_src;

	/* shut up, icc */
	(void)destp;

	predai = (*destap).pred_list.begin();
	predbi = (*destbp).pred_list.begin();
	while(EdgePredPairIterator( destap, destbp, eapp, ebpp, &temp_src))
	{
//		cout << "sourcep="<<sourcep<<" and temp_src="<<temp_src<<endl;
		if (sourcep == temp_src)
			return true;
	}
	cout << " eh\n";
	return false;
}

int CrossProduct::DijkstraBackward (bool global)
{
	
	AlignNode *dest, *src;
	vector<AlignEdge*>::const_iterator ei;
	vector<Node*>::reverse_iterator ai;
	vector<Node*>::reverse_iterator bi;
	int maxscore = MINSCORE;
	//bestend = 0;

	/* shut up, icc */
	(void)global;
	
	for (ai=a.node_list.rbegin(); ai!= a.node_list.rend(); ai++) 
	{
		for (bi=b.node_list.rbegin(); bi!= b.node_list.rend(); bi++)
		{
			src = matrix[(**ai).index][(**bi).index] ;
			for ( ei = (*src).succ_list.begin(); ei != (*src).succ_list.end(); ei++)
			{
				dest = &((**ei).dest);
				if (dest->bscore > MINSCORE && (**ei).score+dest->bscore>src->bscore)
				{
					src->bscore = (**ei).score+dest->bscore;
					//src->dijkstra_parent = dest; 
					if (src->bscore > maxscore)
					{
						maxscore = src->bscore;
					//	bestend = src;
					//	cout << "maxscore=" << maxscore <<" bestend = " << *src << endl;
						
					}
				}
			}
		}
	}
	return  matrix[0][0]->bscore;/// CHANGE THIS	
}


inline bool SmallCrossProduct::EdgePredPairIterator(Node * nap, Node * nbp, Edge ** eapp, Edge ** ebpp, SmallAlignNode **sanpp)
{
	// 	vector<Edge*>::iterator predai, predbi;
	while( predai != ((*nap).pred_list).end()) 
	{ 
		while (predbi != ((*nbp).pred_list).end())
		{
			*eapp = *predai; 
			*ebpp = *predbi;
			*sanpp = &matrix[(**predai).source.index][(**predbi).source.index];
			predbi++;
			return true;
		}
		
		predbi = (*nbp).pred_list.begin();
		predai++;
	}
	return false;
}

inline bool  SmallCrossProduct::VertexPairIterator( Node ** napp,  Node ** nbpp, SmallAlignNode ** sanpp)
{

	/*
	 * vector<Node*>::iterator ai, bi;
	 */ 
	
	while (ai!= a.node_list.end())
	{
		while (bi!= b.node_list.end())
		{
			*napp = *ai; // *ai is Node*, but ai IS NOT Node ** - * is overloaded
			*nbpp = *bi;
			*sanpp = &matrix[(**ai).index][(**bi).index]; //In SmallCrossProduct matrix[i][j] is an element, not a pointer
			bi++;
			return true;
		}
		ai++;
		bi = b.node_list.begin();
	}
	ai = a.node_list.begin();
	return false;
}

SmallCrossProduct::SmallCrossProduct(Graph & a, Graph & b):
		a(a),b(b),column_length(a.column_length + b.column_length)
{
	int asize = a.num_nodes;
	int bsize = b.num_nodes;
	matrix.resize(asize,  vector<SmallAlignNode>(bsize));
	ai = a.node_list.begin();
	bi = b.node_list.begin();
}

void CrossProduct::initscore(bool global)
{
	vector<Node*>::iterator ai;
	vector<Node*>::iterator bi;
	for (ai=a.node_list.begin(); ai!= a.node_list.end(); ai++) 
	{
		for (bi=b.node_list.begin(); bi!= b.node_list.end(); bi++)
		{
			matrix[(**ai).index][(**bi).index]->fscore = (global ? MINSCORE : 0);
			matrix[(**ai).index][(**bi).index]->bscore = (global ? MINSCORE : 0);
			matrix[(**ai).index][(**bi).index]->dijkstra_parent = NULL;
		}
	}
	
	if (global)
	{
		matrix[a.node_list.front()->index][b.node_list.front()->index]->fscore = 0;
		matrix[a.node_list.back()->index][b.node_list.back()->index]->bscore = 0;
	}
		
}
