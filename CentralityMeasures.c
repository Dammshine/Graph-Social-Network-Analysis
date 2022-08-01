// Centrality Measures API implementation
// COMP2521 Assignment 2

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"




// Closeness
double *DijkstraNtimesSol(Graph g);
struct InfoAboutVertexForBC {
	int ver;
	int size;
	int totalPath;
	int containedVerPath;
	double sum;
};
typedef struct InfoAboutVertexForBC *Info;



// Betweeness
double individualBetCen(NodeData **dataAboutGraph, int ver, int size);
void DFSsearch(NodeData *dataAboutDiagram, int i, int j, Info info, int flag);




// Method 1 :
// To solve All-Pair Shortest Paths problem, one can ran single-source shortest path algo for n times
// Complexity would be O(VE * lg(V))
double *DijkstraNtimesSol(Graph g) {
	// Closeness of node x is reciprocal of sum of the lengths 
	// of the shortest paths between node x and other node
	// For more than one connected component, we use Wasserman and Faust formula

	// Creeate array
	int sizeOfVer = GraphNumVertices(g);
	double *clossness = calloc(sizeOfVer, sizeof(double));
	
	// Scan Through array, calculate individual closeness
	for (int i = 0; i < sizeOfVer; i++) {
		// Form single-source shortest path "Diagram"
		NodeData *sp = dijkstra(g, i);

		
		// Calculate closness
		// N = sizeOfVer
		// Sigma all v reachable from u d(u,v)
		int sigma = 0;
		int n = 0;
		for (int i = 0; i < sizeOfVer; i++) {
			if (sp[i].dist != INT_MAX) {
				n++;
				sigma += sp[i].dist;
			}
		}
		if (sigma == 0.0) {
			clossness[i] = 0.0;
			continue;
		}
		double WFfor = ((double) (n - 1) * (n - 1) / (sizeOfVer - 1)) / sigma;
		clossness[i] = WFfor;

		// Free associated resource
		freeNodeData(sp, GraphNumVertices(g));
	}
	return clossness;
}

/*
// Method 2 :
// The Floyd-Warshall algorithm, uses DP formulation
// Complexity would be O(V^3)
double *FWSol(Graph g) {
	return NULL;
}

// Method 3 :
// Johnson's algorithm. This works better for sparse diagram
// Complexity would be O(V^2 *lg(V) + VE)
double *JohnsonSol(Graph g) {
	return NULL;
}
*/


double *closenessCentrality(Graph g) {
	return DijkstraNtimesSol(g);
}

// ===================================================================
//                     Betweenness
// ===================================================================




Info infoCreate(int ver, int size) {
	Info new = malloc(sizeof(*new));
	new->ver = ver;
	new->size = size;
	new->totalPath = 0;
	new->containedVerPath = 0;
	new->sum = 0.0;
	return new;
}

// DFS search
void DFSsearch(NodeData *dataAboutDiagram, int i, int j, Info info, int flag) {
	// Base case
	if (i == j) {
		info->totalPath++;
		// If flag raised, indicate this path passes vertice
		if (flag == 1) {
			info->containedVerPath++;
		}
		return;
	}

	// raise flag
	if (j == info->ver) {
		flag = 1;
	}

	// Scan for viable path from j to i
	struct PredNode *temp = dataAboutDiagram[j].pred;
	while (temp != NULL) {
		DFSsearch(dataAboutDiagram, i, temp->v, info, flag);
		temp = temp->next;
	}
	return;
}

// Calculate Betweeness for a given node
double individualBetCen(NodeData **dataAboutGraph, int ver, int size) {
	Info info = infoCreate(ver, size);

	for (int i = 0; i < size; i++) {
		if (i == ver) continue;
		// Inspect shortest path diagram [diagram i] 
		for (int j = 0; j < size; j++) {
			if (j == ver) continue;
			// Inspect shortest path between i and j in shortest path from [diagram i] 
			// DFS inplementation
			// Skip non-conneted component
			if (dataAboutGraph[i][j].dist == INT_MAX) continue;
			DFSsearch(dataAboutGraph[i], i, j, info, 0);
			// Reset info
			if (info->totalPath != 0) {
				info->sum += (double) info->containedVerPath / info->totalPath;
			}			
			info->totalPath = 0;
			info->containedVerPath = 0;
		}
	}
	double ret = info->sum;
	free(info);
	return ret;
}

// Method 1 :
// Prepare a 2-D array store all shotrest path about each vertex
// Calculate individual betweeness regard to the 2-D array
// Complexity would be O(VE * lg(V))
double *betweennessCentrality(Graph g) {
	double *clossness = calloc(GraphNumVertices(g), sizeof(double));

	// Create 2-D array store all shortest path info about the graph
	NodeData **dataAboutGraph = calloc(GraphNumVertices(g), sizeof(*dataAboutGraph));
	for (int i = 0; i < GraphNumVertices(g); i++) {
		dataAboutGraph[i] = dijkstra(g, i);;
	}

	

	// Calculate individual clossness
	for (int i = 0; i < GraphNumVertices(g); i++) {
		clossness[i] = individualBetCen(dataAboutGraph, i, GraphNumVertices(g));
	}

	for (int i = 0; i < GraphNumVertices(g); i++) {
		freeNodeData(dataAboutGraph[i], GraphNumVertices(g));
	}
	return clossness;
}







struct timeStampVertex {
	int visit;
	int finish;
	int ver;
};
typedef struct timeStampVertex *TimeStampVertex;
bool isDAG(Graph g);
bool isDAGHelper(Graph g, int vertex, int parent, int *visit);
void topologicalSortHelper(Graph g, int ver, TimeStampVertex time, int *t);
int timeCmp(const void *a, const void *b);
int *topologicalSort(Graph g);

// Testing Function (Not part of Assignment)
// ==========================================================
// 					cycle-detection
// ==========================================================
bool isDAG(Graph g) {
	// Using DFS
	int *visit = calloc(GraphNumVertices(g), sizeof(int));
	for (int i = 0; i < GraphNumVertices(g); i++) {
		if (visit[i] == 0) {
			if (isDAGHelper(g, i, i,visit) == false) {
				free(visit);
				return false;
			}
		}
	}
	free(visit);
	return true;
}

bool isDAGHelper(Graph g, int vertex, int parent, int *visit) {
	// MOdify visit status on current vertex
	visit[vertex] = 1;
	AdjList adj = GraphOutIncident(g, vertex);
	AdjList temp = adj;
	while (temp != NULL) {
		// Check status
		if (temp->v == parent) continue;
		else if (visit[temp->v] == 1) {
			return false;
		}
		else {
			if (isDAGHelper(g, temp->v, vertex, visit) == false) {
				return false;
			}
		}
		temp = temp->next;
	}
	return true;
}



// ==========================================================
// 					Topological Sort
// ==========================================================
void topologicalSortHelper(Graph g, int ver, TimeStampVertex time, int *t) {
	// Update visit, update time
	time[ver].visit = *t;
	*t += 1;

	// DFS
	AdjList adj = GraphOutIncident(g, ver);
	AdjList temp = adj;
	while (temp != NULL) {
		// Check status
		if (time[temp->v].visit == 0) {
			topologicalSortHelper(g, temp->v, time, t);
		}
		temp = temp->next;
	}

	// Update finish, update time
	time[ver].finish = *t;
	*t += 1;
}

int timeCmp(const void *a, const void *b) {
	int t1 = ((TimeStampVertex)a)->finish;
	int t2 = ((TimeStampVertex)b)->finish;
	return t2 - t1;
}

int *topologicalSort(Graph g) {
	TimeStampVertex time = calloc(GraphNumVertices(g), sizeof(*time));
	for (int i = 0; i < GraphNumVertices(g); i++) {
		time[i].visit = time[i].finish = 0;
		time[i].ver = i;
	}

	// DFS
	int t = 1;
	for (int i = 0; i < GraphNumVertices(g); i++) {
		if (time[i].visit == 0) {
			topologicalSortHelper(g, i, time, &t);
		}
	}

	// Return the topological order
	int *ret = calloc(GraphNumVertices(g), sizeof(*ret));
	qsort(time, GraphNumVertices(g), sizeof(*time), timeCmp);
	for (int i = 0; i < GraphNumVertices(g); i++) {
		ret[i] = time[i].ver;
	}
	free(time);
	return ret;
}

