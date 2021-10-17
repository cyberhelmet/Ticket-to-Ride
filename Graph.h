#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

typedef struct Edge {
    int to_vertex;
    int weight;
} Edge;

typedef struct EdgeNode {
    Edge edge;
    struct EdgeNode* next;
} EdgeNode;

typedef struct EdgeList {
    EdgeNode* head;
} EdgeList;

typedef struct Graph {
	int vertices;
    EdgeList* edges;
} Graph;

// Function to create graph of give vertices
struct Graph* createGraph(int vertices)
{
	struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
	graph->vertices = vertices;

	// Create an array of adjacency lists.  Size of array will be V
	graph->edges = malloc(vertices * sizeof(struct EdgeList));

	// Initialize each adjacency list as empty by making head as NULL
	for (int i = 0; i < vertices; ++i)
		graph->edges[i].head = NULL;

	return graph;
}

// Function to add edge in graph
void addEdge(Graph* g, int sourceVertex, int destinationVertex, int weight) {

	for (int i = 0; i < g->vertices; i++) {

		if (i == sourceVertex) {

			for (int j = 0; j < g->vertices; j++) {

				if (j == destinationVertex) {

					// Creating new edges
					EdgeNode* newEdge1 = malloc(sizeof(EdgeNode));
					EdgeNode* newEdge2 = malloc(sizeof(EdgeNode));

					// Initializing new edges
					newEdge1->edge.to_vertex = destinationVertex;
					newEdge2->edge.to_vertex = sourceVertex;
					newEdge1->edge.weight = newEdge2->edge.weight = weight;

					// Adding new edges
					newEdge1->next = g->edges[i].head;
					g->edges[i].head = newEdge1;
					newEdge2->next = g->edges[j].head;
					g->edges[j].head = newEdge2;

					printf("\nEdge added successfully!\n");
					return;
				}
			}

			printf("\nDestination vertex not found!\n");
			return;
		}
	}

	printf("\nSource vertex not found!\n");
}

// Function to input graph (By displaying a menu)
Graph* inputGraph() {

	// Taking number of vertices as input
	int vertices;
	printf("\nEnter number of vertices: ");
	scanf_s("%d", &vertices);

	// Creating graph
	Graph* g = createGraph(vertices);

	// Adding edges in graph
	int choice = -1;

	while (choice != 2) {
		printf("\n1 - Add edge.\n2 - Exit.\n\nEnter choice: ");
		scanf_s("%d", &choice);

		while (choice != 1 && choice != 2) {
			printf("\nInvalid choice!\nEnter choice again: ");
			scanf_s("%d", &choice);
		}

		if (choice == 1) {
			int src;
			int dest;
			int weight;

			printf("\nEnter source vertex: ");
			scanf_s("%d", &src);
			printf("Enter destination vertex: ");
			scanf_s("%d", &dest);
			printf("Enter weight: ");
			scanf_s("%d", &weight);

			addEdge(g, src, dest, weight);
		}
	}

	return g;
}

// Function to display graph
void displayGraph(Graph* g) {
	printf("\nDisplaying graph: \n\n");

	for (int i = 0; i < g->vertices; i++) {
		printf("\nVertex-%d:\n", i);

		EdgeNode* currentEdge = g->edges[i].head;
		int j = 0;

		while (currentEdge != NULL) {
			printf("\n\tEdge-%d\n", j);
			printf("\tDestination vertex: %d\n", currentEdge->edge.to_vertex);
			printf("\tWeight: %d\n", currentEdge->edge.weight);
			currentEdge = currentEdge->next;
			j = j + 1;
		}
	}
}

// Structure to represent a min heap node
struct MinHeapNode {
    int v;
    int key;
};

// Structure to represent a min heap
struct MinHeap {
    int size; // Number of heap nodes present currently
    int capacity; // Capacity of min heap
    int* pos; // This is needed for decreaseKey()
    struct MinHeapNode** array;
};

// A utility function to create a new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v, int key) {
    struct MinHeapNode* minHeapNode = (struct MinHeapNode*)malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->key = key;
    return minHeapNode;
}

// A utilit function to create a Min Heap
struct MinHeap* createMinHeap(int capacity) {
    struct MinHeap* minHeap = (struct MinHeap*)malloc(sizeof(struct MinHeap));
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (struct MinHeapNode**)malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b) {
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->key < minHeap->array[smallest]->key)
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->key < minHeap->array[smallest]->key)
        smallest = right;

    if (smallest != idx) {
        // The nodes to be swapped in min heap
        struct MinHeapNode* smallestNode = minHeap->array[smallest];
        struct MinHeapNode* idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is ampty or not
int isEmpty(struct MinHeap* minHeap) {
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap* minHeap) {
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode* root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decrease key value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap* minHeap, int v, int key) {
    // Get the index of v in  heap array
    int i = minHeap->pos[v];

    // Get the node and update its key value
    minHeap->array[i]->key = key;

    // Travel up while the complete tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->key < minHeap->array[(i - 1) / 2]->key) {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap* minHeap, int v) {
    if (minHeap->pos[v] < minHeap->size)
        return true;
    return false;
}

// A utility function used to print the constructed MST
void printMST(int arr[], int key[], int n) {
	int totalCost = 0;
	printf("\nPrim's MST:\n\n");

	for (int i = 1; i < n; ++i) {
		printf("%d - %d\n", arr[i], i);
		totalCost = totalCost + key[i];
	}

	printf("\nTotal cost: %d\n", totalCost);
}

// The main function that constructs Minimum Spanning Tree (By prim's algorithm)
Graph* prims_mst(struct Graph* graph) {
	int vertices = graph->vertices; // Get the number of vertices in graph
	int* parent = malloc(sizeof(int) * vertices); // Array to store constructed MST
	int* key = malloc(sizeof(int) * vertices); // Key values used to pick minimum weight edge in cut

	// minHeap represents set E
	struct MinHeap* minHeap = createMinHeap(vertices);

	// Initialize min heap with all vertices. Key value of
	// all vertices (except 0th vertex) is initially infinite
	for (int v = 1; v < vertices; ++v) {
		parent[v] = -1;
		key[v] = INT_MAX;
		minHeap->array[v] = newMinHeapNode(v, key[v]);
		minHeap->pos[v] = v;
	}

	// Make key value of 0th vertex as 0 so that it
	// is extracted first
	key[0] = 0;
	minHeap->array[0] = newMinHeapNode(0, key[0]);
	minHeap->pos[0] = 0;

	// Initially size of min heap is equal to V
	minHeap->size = vertices;

	// In the following loop, min heap contains all nodes
	// not yet added to MST.
	while (!isEmpty(minHeap)) {
		// Extract the vertex with minimum key value
		struct MinHeapNode* minHeapNode = extractMin(minHeap);
		int u = minHeapNode->v; // Store the extracted vertex number

		// Traverse through all adjacent vertices of u (the extracted
		// vertex) and update their key values
		struct EdgeNode* pCrawl = graph->edges[u].head;
		while (pCrawl != NULL) {
			int v = pCrawl->edge.to_vertex;

			// If v is not yet included in MST and weight of u-v is
			// less than key value of v, then update key value and
			// parent of v
			if (isInMinHeap(minHeap, v) && pCrawl->edge.weight < key[v]) {
				key[v] = pCrawl->edge.weight;
				parent[v] = u;
				decreaseKey(minHeap, v, key[v]);
			}
			pCrawl = pCrawl->next;
		}
	}

	Graph* g = createGraph(vertices);

	for (int i = 1; i < vertices; i++) {
		addEdge(g, parent[i], i, key[i]);
	}

	// Printing MST
	printMST(parent, key, vertices);
	return g;
}

// Some global variables

int* dist;
int* parent;
int* a;	// Stores Set of Vertices in the desired Path from source to destination vertex excluding the latter
int d = 0;	// Number of Vertices in the desired Path from source to destination vertex excluding the latter

// Function to find the vertex with minimum distance value, from the set of vertices not yet included in shortest path tree
int minDistance(int* dist, bool* sptSet, int V) {

	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < V; v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

void find_parent(int* parent, int j) {

	//Base Case : If j is source
	if (parent[j] == -1)
		return;

	find_parent(parent, parent[j]);
	a[d] = j;
	d++;
}

// Function to print shortest path from source to j using parent array
void printPath(int* parent, int j) {

	if (parent[j] == -1) {
		return;
	}
	printPath(parent, parent[j]);
	printf("%d ", j);
}

// A utility function to print the constructed distance array
int printSolution(int* dist, int n, int x, int V) {
	int src = x;
	printf("Vertex \tDistance from Source\t Path");

	for (int i = 1; i < V; i++) {
		printf("\n%d -> %d \t\t %d\t\t%d ",
			src, i, dist[i], src);
		printPath(parent, i);
	}
	printf("\n");
}

// Function that implements Dijkstra's single source shortest path algorithm for a graph
void dijkstra(int** graph, int src, int V) {

	// The output array. dist[i] will hold the shortest distance from src to i 
	bool* sptSet = (bool*)malloc(sizeof(bool) * V); // sptSet[i] will be true if vertex i is included in 
	// shortest  path tree or shortest distance from src to i is finalized 

	// Initialize all distances as INFINITE and stpSet[] as false 
	for (int i = 0; i < V; i++)
		dist[i] = INT_MAX, sptSet[i] = false, parent[i] = -1;

	// Distance of source vertex from itself is always 0 
	dist[src] = 0;

	// Find shortest path for all vertices 
	for (int count = 0; count < V; count++) {

		// Pick minimum distance vertex from the set of vertices not yet processed.
		// u is always equal to src in the first iteration. 
		int u = minDistance(dist, sptSet, V);

		// Mark the picked vertex as processed 
		sptSet[u] = true;

		// Update dist value of the adjacent vertices of the picked vertex. 
		for (int v = 0; v <= V; v++)

			// Update dist[v] only if is not in sptSet, there is an edge from 
			// u to v, and total weight of path from src to v through u is 
			// smaller than current value of dist[v] 
			if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
				&& dist[u] + graph[u][v] < dist[v]) {
				parent[v] = u;
				dist[v] = dist[u] + graph[u][v];
			}
	}
}

// Striener tree implementation (Using Dijkstra's single source shortest path algorithm)
// Takes array of ticket as parameter and returns MST with minimum cost
void strienerTree(Graph* g, int tickets[], int K) {
	int V = g->vertices;

	int* T = malloc(sizeof(int) * V);	// Set of Vertices in MST T
	int y = 0;	// Number of Vertices in MST T

	dist = malloc(sizeof(int) * V);
	parent = malloc(sizeof(int) * V);
	a = malloc(sizeof(int) * V);

	// Initilizing T and a arrays
	for (int i = 0; i < V; i++) { a[i] = 0; T[i] = 0; }

	int i = 0, j = 0;
	int* remaining_vertices = malloc(sizeof(int) * K);

	// Initializing remaining_vertices array
	for (int i = 0; i < K; i++) { remaining_vertices[i] = 0; }

	int source_vertex = 0;
	int* is_processed = malloc(sizeof(int) * V);

	// Initializing is_processed array
	for (int i = 0; i < V; i++) { is_processed[i] = 0; }

	int** is_added = (int**)malloc(sizeof(int*) * V);
	int** is_included = (int**)malloc(sizeof(int*) * V);

	for (int i = 0; i < V; i++) {
		is_added[i] = (int*)malloc(sizeof(int) * V);
		is_included[i] = (int*)malloc(sizeof(int) * V);
	}

	// Initializing is_added and is_included arrays
	for (int i = 0; i < V; i++) {

		for (int j = 0; j < V; j++) {
			is_added[i][j] = 0;
			is_included[i][j] = 0;
		}
	}

	int total_cost = 0, count = 0, cost = 0, x = 0;

	// Makng a graph array (Matrix)
	int** graph = (int**)malloc(sizeof(int*) * V);

	for (int i = 0; i < V; i++) {
		graph[i] = (int*)malloc(sizeof(int) * V);
	}

	// Initializing graph array
	for (int i = 0; i < V; i++) {

		for (int j = 0; j < V; j++) {
			graph[i][j] = 0;
		}
	}

	// Filling graph matrix
	for (int i = 0; i < V; i++) {
		EdgeNode* currentEdge = g->edges[i].head;

		while (currentEdge != NULL) {
			graph[i][currentEdge->edge.to_vertex] = currentEdge->edge.weight;
			currentEdge = currentEdge->next;
		}
	}

	/*
	 * Step 1 of the algorithm
	 * First Terminal Vertex is added to T
	 * "Start with a subtree T consisting of
	 * one given terminal vertex"
	 */

	dijkstra(graph, tickets[0], V);
	is_processed[tickets[0]] = 1;
	T[0] = tickets[0];
	y++;

	/*
	 * Step 2 of the algorithm starts here
	 * "While T does not span all terminals"
	 */
	count = 1;
	while (count < K)
	{
		/*
		* Step 2 a) of the algorithm
		* "Select a terminal x not in T that is closest
		* to a vertex in T"
		*/

		x = 0;				/* x --> Next Terminal Vertex */
		int min = INT_MAX;

		for (i = 1; i < K; i++)
		{
			if ((min > dist[tickets[i]]) &&
				(dist[tickets[i]] != 0) &&
				((is_processed[tickets[i]]) == 0))
			{
				min = dist[tickets[i]];
				x = tickets[i];
			}
		}

		/*
		* Step 2 b) of the algorithm starts here
		* "Finding Vertex in T which is closest to
		* Next Terminal Vertex to be added"
		*/

		int min_cost = INT_MAX;
		for (i = 0; i < y; i++) {
			cost = 0;
			dijkstra(graph, T[i], V);
			d = 0;
			find_parent(parent, x);

			for (j = 0; j < d; j++) {

				if (j == 0) {

					if (is_added[T[i]][a[0]] == 0) {
						is_added[T[i]][a[0]] = 1;
						cost += graph[T[i]][a[0]];
					}
				}

				else {

					if (is_added[a[j - 1]][a[j]] == 0) {
						is_added[a[j - 1]][a[j]] = 1;
						cost += graph[a[j - 1]][a[j]];
					}
				}
			}

			if (cost < min_cost) {
				min_cost = cost;
				source_vertex = T[i];
			}

			for (int l = 0; l < V; l++) {

				for (int m = 0; m < V; m++) {
					is_added[l][m] = is_included[l][m];
				}
			}
		}

		/*
		* x is connected with "source_vertex" in T
		* "Adding to T the shortest path that connects x with T"
		* Total Cost is calculated here
		*/

		dijkstra(graph, source_vertex, V);

		d = 0;
		find_parent(parent, x);

		for (j = 0; j < d; j++) {

			if (j == 0) {

				if (is_included[source_vertex][a[0]] == 0) {
					is_included[source_vertex][a[0]] = 1;
					total_cost += graph[source_vertex][a[0]];

					if (is_processed[a[0]] == 0) {
						is_processed[a[0]] = 1;
						T[y] = a[0];
						y++;
					}
				}
			}

			else {

				if (is_included[a[j - 1]][a[j]] == 0) {
					is_included[a[j - 1]][a[j]] = 1;
					total_cost += graph[a[j - 1]][a[j]];

					if (is_processed[a[j - 1]] == 0) {
						is_processed[a[j - 1]] = 1;
						T[y] = a[j - 1];
						y++;
					}

					if (is_processed[a[j]] == 0) {
						is_processed[a[j]] = 1;
						T[y] = a[j];
						y++;
					}
				}
			}
		}
		count++;
	}

	// Displaying all edges in minimum cost MST
	printf("\nIn striener tree function:\n\n");

	for (i = 0; i < V; i++) {

		for (j = 0; j < V; j++) {

			if (is_included[i][j] == 1) {
				printf("%d - %d\n", i, j);
			}
		}
	}

	printf("\nTotal Cost : %d\n", total_cost);
}