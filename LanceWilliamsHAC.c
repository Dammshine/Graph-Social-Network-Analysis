// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY 	DBL_MAX
#define ROOT_FLAG 	-1


double minDouble(double a, double b);
double maxDouble(double a, double b);
int maxInt(int a, int b);
int minInt(int a, int b);


// Return max, min among two int
int minInt(int a, int b) {
	if (a < b) return a;
	return b;
}

int maxInt(int a, int b) {
	if (b <= a) return a;
	return b;
}


// Return max, min among
// If one of a,b is inf, return valid one
double maxDouble(double a, double b) {
	if (a == INFINITY && b == INFINITY) return INFINITY;
	if (a == INFINITY) return b;
	if (b == INFINITY) return a;
	return (a > b) ? a : b;
}

double minDouble(double a, double b) {
	if (a == INFINITY && b == INFINITY) return INFINITY;
	if (a == INFINITY) return b;
	if (b == INFINITY) return a;
	return (a > b) ? b : a;
}

// Merge two index in Dendrogram
void DendrogramMerge(Dendrogram *dtree, int a, int b, int size) {
	// Sentinel for valid index
	if (a >= size || b >= size) return;
	if (a < 0 || b < 0) return;


	// Always move the smaller to left, move greater one on right
	// If merging empty node, put empty one on the right
	int min = minInt(a, b);
	int max = maxInt(a, b);
	Dendrogram parentNode = malloc(sizeof(*parentNode));
	parentNode->vertex = INT_MAX;

	// Using -1 as a flag, indicate returned Dendrogram
	if (size == 2) parentNode->vertex = ROOT_FLAG;
	parentNode->left = dtree[min];
	parentNode->right = dtree[max];

	// Replace both link to the common acestor
	dtree[min] = parentNode;
	for (int i = max + 1; i < size; i++) {
		dtree[i - 1] = dtree[i]; 
	}
}

// =====================================================================
// 							Matrix
// =====================================================================
void matrixMerge(double **curr, double **prev, Dendrogram *denTree, int size, int mode);
void matrixFree(double **matrix, int size);
double **matrixCreate(int size);
void initialize(double **matrix, int size, Graph g);
void matrixCopy(double **matrixA, double **matrixB, int size);


// Struct used for finding merged col and row
struct coordinate {
	int i;
	int j;
	double weight;
};

// The Dendrogram will denote the index for the curr matirx
// prev matrix used in the process of merge
void matrixMerge(double **curr, double **prev, Dendrogram *denTree, int size, int mode) {
	// Use coordinate to choose the smallest from the graph
	struct coordinate coor;
	coor.i = -1;
	coor.j = -1;
	coor.weight = DBL_MAX;

	// Scan to find the smallest which is also the most important link
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (curr[i][j] < coor.weight) {
				coor.weight = curr[i][j];
				coor.i = i;
				coor.j = j;
			} 
		}
	}

	// Now merge denTree
	DendrogramMerge(denTree, coor.i, coor.j, size);
	
	// Saving curr into prev
	matrixCopy(prev, curr, size);

	// min, and max denote the smaller/greater cell in merged cluster
	// new matrix will merge two cluster into min cluster
	int min = minInt(coor.i, coor.j);
	int max = maxInt(coor.i, coor.j);

	// Update curr matrix
	// Map cell in previous matrix to the new matrix
	int row = 0;
	int col = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// Case for non-valid cell
			if (i == max && j == max) curr[min][min] = INFINITY;
			else if (i == j) curr[row][col] = INFINITY;

			// Case for general cell (Do not associate with updated cell)
			if ((i != max && i != min) && (j != max && j != min)) {
				curr[row][col] = prev[i][j];
			} // Case for merged cell (associated with min and max)
			else if (i == max || j == max) { // Case 1 : If cell associated with max, ignore
				;
			} // Case 2: Update cell associated with min cell
			else {
				// If col is min, replace it with mode specific 
				if (i == min) {
					if (mode == SINGLE_LINKAGE) curr[row][col] = minDouble(prev[min][j], prev[max][j]);
					else curr[row][col] = maxDouble(prev[min][j], prev[max][j]);
				} // If row is min
				else {
					if (mode == SINGLE_LINKAGE) curr[row][col] = minDouble(prev[i][min], prev[i][max]);
					else curr[row][col] = maxDouble(prev[i][min], prev[i][max]);
				}
			}

			// For ignore the max in matrix
			if (j != max) col++;
		}
		col = 0;
		if (i != max) row++;
	}
}

// Free associated memory of matrix
void matrixFree(double **matrix, int size) {
	for (int i = 0; i < size; i++) {
		free(matrix[i]);
	}
	free(matrix);
}

// Create matrix, initialize all grid as INFINITY
double **matrixCreate(int size) {
	double **matrix = calloc(size, sizeof(*matrix));
	for (int i = 0; i < size ; i++) {
		matrix[i] = calloc(size, sizeof(**matrix));
		for (int j = 0; j < size ; j++)  {
			matrix[i][j] = INFINITY;
		}
	}
	return matrix;
}


// Initialize matrix's distance, as represended by the graph g
void initialize(double **matrix, int size, Graph g) {
	for (int i = 0; i < size; i++) {
		// First initiaze all node that have out edge
		// Using new dist conversion 1 / weight
		AdjList temp = GraphOutIncident(g, i);
		while (temp != NULL) {
			matrix[i][temp->v] = (double) 1 / temp->weight;
			temp = temp->next;
		}

		// Then add new, or relax previous outEdges
		temp = GraphInIncident(g, i);
		while (temp != NULL) {
			double inEdgeWeight = (double) 1 / temp->weight;
			if (inEdgeWeight < matrix[i][temp->v]) {
				matrix[i][temp->v] = inEdgeWeight;
			}
			temp = temp->next;
		}
	}
}

// Copy matrixB to matrixA, for size * size
void matrixCopy(double **matrixA, double **matrixB, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			matrixA[i][j] = matrixB[i][j];
		}
	}
}

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
	// Initialize denGrogram
	Dendrogram *denTree = calloc(GraphNumVertices(g), sizeof(*denTree));
	
	// Create Demdrogram for each of vertex
	int size = GraphNumVertices(g);
	for (int i = 0; i < GraphNumVertices(g); i++) {
		denTree[i] = malloc(sizeof(*denTree[i]));
		denTree[i]->left = NULL;
		denTree[i]->right = NULL;
		denTree[i]->vertex = i;
	}

	// Create a matrix I guess
	double **currMatrix = matrixCreate(GraphNumVertices(g));
	double **prevMatrix = matrixCreate(GraphNumVertices(g));

	// Initialize curr Matrix
	initialize(currMatrix, GraphNumVertices(g),  g);

	// While size allowed, merge two cluster in currMatrix
	while (size > 0) {
		matrixMerge(currMatrix, prevMatrix, denTree, size, method);
		size--;
	}

	// Retrive the root of Dendrogram
	Dendrogram ret = NULL;
	for (int i = 0; i < GraphNumVertices(g); i++) {
		if (denTree[i]->vertex == ROOT_FLAG) ret = denTree[i];
	}

	// Free associate memory
	matrixFree(prevMatrix, GraphNumVertices(g));
	matrixFree(currMatrix, GraphNumVertices(g));
	free(denTree);

	return ret;
}

