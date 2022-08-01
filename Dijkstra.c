// Dijkstra API implementation
// COMP2521 Assignment 2

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"


// ========================================================
// 			Helper function for PredNode
// ========================================================
// Create PredNode contain member v and link
PredNode *PredNodeCreate(Vertex v, PredNode *link);


// Insert Vertex into PredNode linked list in acending order
PredNode *PredNodeInsert(PredNode *n, Vertex v);


PredNode *PredNodeEmplace(PredNode *n, Vertex v);

void PredNodeFree(PredNode *n);



// ========================================================
// 			Helper function for NodeDate
// ========================================================
// Garunteed to inset Vertex v in correct position inside of n's pred
void NodeDataInsert(NodeData n, Vertex v) {
	n.pred = PredNodeInsert(n.pred, v);
}


// Also it assumed the edge weight are positive, however don't get mentioned in graph ADT
/**
 * Finds  the shortest paths from a given source node to all other nodes
 * in the given graph. See the spec for more details.
 *
 * The  function  returns  a  NodeData  array,  with length equal to the
 * number of nodes in the graph, where index v of the array contains the
 * distance and predecessor list for node v. The predecessor lists  must
 * be in ascending order of vertex number.
 */

// Insert reachable/valid node from v into PQ
// Update the distance if need
void AdjenctVerticeInsert(Graph g, Vertex v, PQ pq, NodeData *a, char *visit) {
	// Inspect curr node
	int currWeight = a[v].dist;

	// Get outEdge collection
	AdjList outedges = GraphOutIncident(g, v);
	AdjList temp = outedges;
	//printf("Currenly visit %d ", v);
	while (temp != NULL) {
		//printf("Inspect Edge %d with %c ", temp->v, visit[temp->v]);
		// Scan through viable edge
		int inspect = temp->v;
		int inspectWeight = temp->weight;
		inspectWeight += currWeight;

		// Skip visited node
		if (visit[inspect] == 'b') {
			temp = temp->next;
			continue;
		}

		// Update inspected edge, if it's unreachable before
		if (visit[inspect] == 'w') {
			a[inspect].dist = inspectWeight;
			a[inspect].pred = PredNodeInsert(a[inspect].pred, v);
			visit[inspect] = 'g';
			PQInsert(pq, inspect, inspectWeight);
		}

		// Update the distance for reachable vertex
		else if (inspectWeight < a[inspect].dist) {
			a[inspect].pred = PredNodeEmplace(a[inspect].pred, v);
			a[inspect].dist = inspectWeight;
			PQUpdate(pq, inspect, inspectWeight);
		}
		else if (inspectWeight == a[inspect].dist) {
			a[inspect].pred = PredNodeInsert(a[inspect].pred, v);
		}
		temp = temp->next;
	}
	//printf("\n");
	
	// Blacken out current vertex
	visit[v] = 'b';
}


NodeData *dijkstra(Graph g, Vertex src) {
	// Create NodeData array
	int graphSize = GraphNumVertices(g);
	NodeData *nodeDataArray = calloc(graphSize, sizeof(NodeData));
	for (int i = 0; i < graphSize; i++) {
		nodeDataArray[i].dist = INT_MAX;
		nodeDataArray[i].pred = NULL;
	}
	nodeDataArray[src].dist = 0;
	nodeDataArray[src].pred = NULL;

	// Create array to track visit status
	// 'w' white : unreachable
	// 'g' grey  : reachable by previous visited node
	// 'b' black : visited by dijkstra
	char *visit = calloc(graphSize, sizeof(char));
	for (int i = 0; i < graphSize; i++) {
		visit[i] = 'w';
	}
	visit[src] = 'g';

	// Create priority Queue to track all current distance
	PQ pq = PQNew();
	PQInsert(pq, src, 0);

	while (PQIsEmpty(pq) != true) {
		Vertex inspect = PQDequeue(pq);
		AdjenctVerticeInsert(g, inspect, pq, nodeDataArray, visit);
	}

	PQFree(pq);
	free(visit);
	return nodeDataArray;
}



void freeNodeData(NodeData *data, int nV) {
	for (int i = 0; i < nV; i++) {
		PredNodeFree(data[i].pred);
	}
	free(data);
}



// Insert PredNode into correct position
PredNode *PredNodeInsert(PredNode *n, Vertex v) {
	// Base case, insert at the end
	if (n == NULL) return PredNodeCreate(v, NULL);

	// Base case 2: v < current node, v need to be insert before current node
	if (v < n->v) {
		PredNode *temp = PredNodeCreate(v, n);
		return temp;
	}	// Recursive case: advance current node
	else {
		n->next = PredNodeInsert(n->next, v);
		return n;
	}
}

// Create PredNode 
PredNode *PredNodeCreate(Vertex v, PredNode *link) {
	PredNode *pn = malloc(sizeof(*pn));
	pn->v = v;
	pn->next = link;
	return pn;
}

void PredNodeFree(PredNode *n) {
	while (n != NULL) {
		PredNode *temp = n;
		n = n->next;
		free(temp);
	}
}

// Replace PredNode with new Node
PredNode *PredNodeEmplace(PredNode *n, Vertex v) {
	PredNodeFree(n);
	return PredNodeCreate(v, NULL);
}
