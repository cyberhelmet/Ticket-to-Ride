#include "Graph.h"

int main() {

	// Taking graph as input
	Graph* g = inputGraph();

	// Displaying graph
	displayGraph(g);

	// Prim's algorithm
	prims_mst(g);

	// Taking number of cities as input
	int cities;
	printf("\nEnter number of cities for ticket: ");
	scanf_s("%d", &cities);

	int* ticket = (int*)malloc(sizeof(int) * cities);

	for (int i = 0; i < cities; i++) {
		printf("\nEnter city-%d: ", i + 1);
		scanf_s("%d", &ticket[i]);
	}

	// Calling striener tree function
	strienerTree(g, ticket, cities);

	return 0;
}