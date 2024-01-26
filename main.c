#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


// This program finds communities within a graph utilizin edge betweenness values.
// Graph is represented by adjacency lists and graph type is undirected.
// Edge betweenness degree of an edge is how much a shortest path from arbitrary 2 nodes passes from this edge. 
// BFS is used for finding all possible shortest paths from a node to another nodes.
// BFS is implemented using a queue. After finding all shortest paths each edge value forming this path is increased by 1.
// After updating edge betweenness applying BFS to each node in graph then we find the highest degree edge.
// We remove this edge from graph and detect communities by finding connected components using DFS.
// We keep iterating above steps to loosen the edge connections in the graph
// Until a predefined minimum community size has been reached or iterations has been exceeded.



// Necessary queue structure and functions for implementing BFS.
struct Queue {
    int front, rear, size;
    unsigned capacity;
    int* array;
};

struct Queue* createQueue(unsigned capacity) {
    struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));
    queue->capacity = capacity;
    queue->front = queue->size = 0;
    queue->rear = capacity - 1;
    queue->array = (int*)malloc(queue->capacity * sizeof(int));
    return queue;
}

int isFull(struct Queue* queue) {
    return (queue->size == queue->capacity);
}

bool isEmpty(struct Queue* queue) {
    return (queue->size == 0);
}

void enqueue(struct Queue* queue, int item) {
    if (isFull(queue))
        return;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->array[queue->rear] = item;
    queue->size = queue->size + 1;
}

int dequeue(struct Queue* queue) {
    if (isEmpty(queue))
        return INT_MIN;
    int item = queue->array[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}

// A structure to represent an adjacency list node
struct AdjListNode {
    int dest;
    struct AdjListNode* next;
};
 
// A structure to represent an adjacency liat
struct AdjList {
    struct AdjListNode *head; // pointer to head node of list
};
 
// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V and taken as input (number of vertices in graph)
struct Graph {
    int V;
    struct AdjList* array;
};
 
// A utility function to create a new adjacency list node
struct AdjListNode* newAdjListNode(int dest) {
    struct AdjListNode* newNode = (struct AdjListNode*) malloc(
            sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}
 
// A utility function that creates a graph of V vertices
struct Graph* createGraph(int V) {
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;
 
    // Create an array of adjacency lists.  Size of array will be V
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));
 
    // Initialize each adjacency list as empty without edges.
    int i;
    for (i = 0; i < V; ++i)
        graph->array[i].head = NULL;
 
    return graph;
}
 
// Function to add an edge given two vertices (A,B) as parameters.
void addEdge(struct Graph* graph, int src, int dest) {
    // Add an edge from src to dest.  A new node is added to the adjacency
    struct AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
 
    // Since graph is undirected, add an edge from dest to src also
    newNode = newAdjListNode(src);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
}


// Function to remove edge function takes an edge v1,v2 (src, dest).
// Since graph is undirected we also remove v2,v1.
void removeEdge(struct Graph* graph, int src, int dest) {
    // Remove edge from src to dest
    struct AdjListNode *temp = graph->array[src].head, *prev = NULL;
    while (temp != NULL) {
        if (temp->dest == dest) {
            if (prev == NULL) {
                graph->array[src].head = temp->next;
            } else {
                prev->next = temp->next;
            }
            struct AdjListNode* toFree = temp;
            temp = temp->next; // Move to next node before freeing the current node
            free(toFree);
            continue;
        }
        prev = temp;
        temp = temp->next;
    }

    // Remove edge from dest to src
    temp = graph->array[dest].head;
    prev = NULL;
    while (temp != NULL) {
        if (temp->dest == src) {
            if (prev == NULL) {
                graph->array[dest].head = temp->next;
            } else {
                prev->next = temp->next;
            }
            struct AdjListNode* toFree = temp;
            temp = temp->next;
            free(toFree);
            continue;
        }
        prev = temp;
        temp = temp->next;
    }
}
 
// A utility function to print the adjacenncy list representation of graph
void printGraph(struct Graph* graph) {
    int v;
    for (v = 0; v < graph->V; ++v) {
        struct AdjListNode* pCrawl = graph->array[v].head;
        printf("\n Adjacency list of vertex %d\n head ", v);
        while (pCrawl) {
            printf("-> %d", pCrawl->dest);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}
 
 
// This function performs BFS and calculates all shortest paths from a node A to all other nodes.
// distance array represents the length of the shortest path from startVertex to another node v2 distance[v2].
// Initially all distances to another nodes are infinite and distance[startVertex] = 0.
// We use predecessors matrix which holds adjacent nodes to a node that is in the shortest path.
// We can select all adjacent nodes branching from a vertex for multiple paths as long as minimum distance is preserved.
// Later on we backtrack from predecessors matrix to find all possible shortest paths that is leading to startVertex. 
int** BFS_AllShortestPaths(struct Graph* graph, int startVertex) {
    int i, j;
    int* distance = (int*)malloc(graph->V * sizeof(int));
    int** predecessors = (int**)malloc(graph->V * sizeof(int*));
    for (i = 0; i < graph->V; i++) {
        distance[i] = INT_MAX;
        predecessors[i] = (int*)malloc(graph->V * sizeof(int));
        for (j = 0; j < graph->V; j++) {
            predecessors[i][j] = -1;
        }
    }

    struct Queue* queue = createQueue(graph->V);
    distance[startVertex] = 0;
    enqueue(queue, startVertex);
    //printf("Starting BFS from Vertex %d\n", startVertex);

    while (!isEmpty(queue)) {
        int currentVertex = dequeue(queue);
        //printf("Dequeued Vertex %d\n", currentVertex);

        struct AdjListNode* temp = graph->array[currentVertex].head;
        while (temp) {
            int adjVertex = temp->dest;
            //printf("Checking Vertex %d from Vertex %d\n", adjVertex, currentVertex);

            if (distance[adjVertex] > distance[currentVertex] + 1) {
                distance[adjVertex] = distance[currentVertex] + 1;
                predecessors[adjVertex][0] = currentVertex;
                enqueue(queue, adjVertex);
                //printf("Updated Distance for Vertex %d is now %d\n", adjVertex, distance[adjVertex]);
                //printf("Updated Predecessors for Vertex %d: ", adjVertex);
                for (j = 0; j < graph->V && predecessors[adjVertex][j] != -1; j++) {
                    //printf("[%d] %d ", j, predecessors[adjVertex][j]);
                }
                //printf("\n");
            } else if (distance[adjVertex] == distance[currentVertex] + 1) {
                int index = 0;
                while (predecessors[adjVertex][index] != -1) index++;
                predecessors[adjVertex][index] = currentVertex;
                //printf("Added additional predecessor %d for Vertex %d at index %d\n", currentVertex, adjVertex, index);
            }

            temp = temp->next;
        }
    }
    
	/*
    printf("BFS Complete. Distance and Predecessors matrix:\n");
    
    printf("\nDistances array:\n");
    for(i=0;i<graph->V;i++){
    	printf("%d ",distance[i]);
	}
	
    printf("\n\nPredecessors Matrix:");
    
    for(i=0;i<graph->V;i++){
    	printf("\n");
    	for(j=0;j<graph->V;j++){
    		printf("%d ",predecessors[i][j]);
		}
	}
    
    for (i = 0; i < graph->V; i++) {
        printf("\nVertex %d: Distance: %d, Predecessors: ", i, distance[i]);
        for (j = 0; j < graph->V && predecessors[i][j] != -1; j++) {
            printf("[%d] %d ", j, predecessors[i][j]);
        }
        printf("\n");
    }
    
    */
	return predecessors;
}


// Takes an edgeBetweenness matrix and updates it for a vertex by backtracking predecessors matrix derived from BFS.
// When an edge (v1,v2) is passed while backtracking then edgeBetweenness[v1][v2]++, edgeBetweenness[v2][v1]++
void updateEdgeBetweenness(struct Graph* graph, int startVertex, int** predecessors, int** edgeBetweenness) {
    int i, endVertex;
    
    //printf("Start vertex: %d\n",startVertex);

    for (endVertex = 0; endVertex < graph->V; endVertex++) {
        if (endVertex != startVertex) {
            //printf("%d. nodedan %d. nodea giden yollar:\n", startVertex, endVertex);

            for (i = 0; predecessors[endVertex][i] != -1; i++) {
                int pred = predecessors[endVertex][i];
                int currentVertex = endVertex;

                //printf("%d", endVertex);

                while (currentVertex != startVertex) {
                    //printf(" -> %d", pred);

                    // Increment edge betweenness for the edge from pred to currentVertex
                    edgeBetweenness[pred][currentVertex] += 1;
                    edgeBetweenness[currentVertex][pred] += 1; // For undirected graphs, remove for directed graphs

                    // Move to the next vertex in the path
                    currentVertex = pred;
                    if (currentVertex != startVertex) {
                        pred = predecessors[pred][0]; // Using the first predecessor for simplicity
                    }
                }

                //printf("\n");
            }
        }
    }
    //printf("\n------------\n");
}


// Performs updateEdgeBetweenness function for all nodes in the graph and gets the edgeBetweenness matrix.
// Highest degree edge can be found and removed later on using removeEdge function.
int** getEdgeBetweenness(struct Graph* graph){
	int i,j;
	int** edgeBetweenness = (int**)malloc(graph->V * sizeof(int*));
	for (i = 0; i < graph->V; i++) {
        edgeBetweenness[i] = (int*)malloc(graph->V * sizeof(int));
        for (j = 0; j < graph->V; j++) {
            edgeBetweenness[i][j] = 0;
        }
    }
    
    for(i=0; i<graph->V;i++){
    	
    	int startVertex = i;
    	int** predecessors = BFS_AllShortestPaths(graph, startVertex);
    	updateEdgeBetweenness(graph, startVertex, predecessors, edgeBetweenness);
    	
	}
	
	
	return edgeBetweenness;
    
}


// Utilizing DFS for finding connected components (communities)

void DFS(struct Graph* graph, int v, bool visited[], int* community, int* size){
	visited[v] = true;
	community[*size] = v;
	(*size)++;
	
	struct AdjListNode* temp = graph->array[v].head;
	while (temp) {
		int adjVertex = temp->dest;
		if(!visited[adjVertex]){
			DFS(graph, adjVertex, visited, community, size);
		}
		temp = temp->next;
	}
}

struct CommunityInfo {
    int minCommunitySize;
    int communityCount;
}CommunityInfo;


// maxIterations = threshold for iterations without community count change.
// minSize = When a community size gets below minSize we stop iterating.


struct CommunityInfo detectCommunities(struct Graph* graph, int maxIterations, int minSize){
	int i, j, v;
	int iterations = 0;
	struct CommunityInfo communityInfo;
	communityInfo.minCommunitySize = graph->V; // Max possible size
    communityInfo.communityCount = 0;
    int prevCommunityCount = 0;
    
	
	while(communityInfo.minCommunitySize > minSize && iterations < maxIterations){
		prevCommunityCount = communityInfo.communityCount;
		communityInfo.communityCount = 0;
		printf("\n============================================\n");
		
		printf("\nIteration: %d\n\n",iterations+1);
		
		bool* visited = (bool*)malloc(graph->V * sizeof(bool));
		for(i=0; i<graph->V; i++){
			visited[i] = false;
		}
		
		// First find communities, community count, community size using DFS.
		int* community = (int*)malloc(graph->V * sizeof(int));
		int size;
		
		for(v=0; v<graph->V; v++){
			
			if(!visited[v]){
				size = 0;
				DFS(graph, v, visited, community, &size);
				printf("\nCommunity %d: ",communityInfo.communityCount);
				printCommunity(community, size);
				communityInfo.communityCount++;
				if(size < communityInfo.minCommunitySize){
					communityInfo.minCommunitySize = size;
				}
			}
		}
        
        // When there is no change in community count increment iterations else reset it.
        if(communityInfo.communityCount == prevCommunityCount){
			iterations++;
		} else{
			iterations = 0;
		}
		
        printf("\n============================================\n");
		free(visited);
		free(community);
		
		if(communityInfo.minCommunitySize > minSize && iterations < maxIterations){
			// Community conditions are not satisfied so we remove the edge with highest degree.
			int** edgeBetweenness = getEdgeBetweenness(graph);
			int highestBetweenness = 0;
			int v1 = -1, v2 = -1;
			for (i = 0; i < graph->V; i++) {
	            for (j = i + 1; j < graph->V; j++) { // Assuming undirected graph
	                if (edgeBetweenness[i][j] > highestBetweenness) {
	                    highestBetweenness = edgeBetweenness[i][j];
	                    v1 = i;
	                    v2 = j;
	                }
	            }
	        }
	        printf("\nEdge betweenness matrix: \n");
	        for (i=0; i<graph->V; i++){
	        	printf("\n");
	        	for(j=0; j<graph->V; j++){
	        		printf("%d ",edgeBetweenness[i][j]);
				}
			}
	        // Remove the edge with the highest betweenness and keep iterating.
	        if (v1 != -1 && v2 != -1) {
	            removeEdge(graph, v1, v2);
	            communityInfo.minCommunitySize = graph->V;
	            printf("\n\nRemoved edge with highest betweenness: (%d, %d)\n", v1, v2);
	        }
	        
		}
		
	}
	
	if(communityInfo.minCommunitySize <= minSize){
		printf("Reached minimum community size.\n");
	} else {
		printf("Reached maximum iterations.\n");
	}
	
	return communityInfo;
	
}



// Function to print the community
void printCommunity(int community[], int size) {
	int i;
    for (i = 0; i < size; i++) {
        printf("%d ", community[i]);
    }
    printf("\n");
}



int main() {
	printf("Please enter input file txt path: ");
	char filename[100];
    gets(filename);
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("\nError opening file");
        return -1;
    }

    int V;
    // Read graph from an input txt file.
    // Format: 
    // 8; (first line vertex count)
    // 0:1, (Adjacency list of vertex 0)
    // 1:0,2; (Adjacency list of vertex 1)
    // 2:1; (Adjacency list of vertex 2)
    // ..
    // ..
    // 7:4,6;
    printf("\n");
    fscanf(file, "%d;", &V);  // Read vertex count
    struct Graph* graph = createGraph(V);

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        int src, dest;
        char* token = strtok(line, ":,");
        src = atoi(token);

        while ((token = strtok(NULL, ":,")) != NULL) {
            dest = atoi(token);
            addEdge(graph, src, dest);
        }
    }
    fclose(file);
    
    
    printGraph(graph);
    int maxIterations, minCommSize;
    printf("\nPlease enter max iterations: ");
    scanf("%d", &maxIterations);
    printf("\nPlease enter min community size: ");
    scanf("%d", &minCommSize);
    
    
    
    
    
    


	struct CommunityInfo communityInfo;
	communityInfo = detectCommunities(graph,maxIterations, minCommSize);
	
	printf("\nmin community size: %d\n", communityInfo.minCommunitySize);
	
	printf("community count: %d\n", communityInfo.communityCount);
	

	
	/*
	int** predecessors = BFS_AllShortestPaths(graph, 0);
	
	int** edgeBetweenness = (int**)malloc(graph->V * sizeof(int*));
	for (i = 0; i < graph->V; i++) {
        edgeBetweenness[i] = (int*)malloc(graph->V * sizeof(int));
        for (j = 0; j < graph->V; j++) {
            edgeBetweenness[i][j] = 0;
        }
    }
	
	updateEdgeBetweenness(graph, 0, predecessors, edgeBetweenness);
	
	for (i=0; i<graph->V; i++){
		printf("\n");
		for(j=0;j<graph->V;j++){
			printf("%d ",edgeBetweenness[i][j]);
		}
	}
	*/
	
	/*
    V = 8;
    struct Graph* graph = createGraph(V);
    addEdge(graph, 0, 1);
    addEdge(graph, 1, 2);
    addEdge(graph, 2, 3);
    addEdge(graph, 3, 4);
    addEdge(graph, 4, 5);
    addEdge(graph, 2, 6);
    addEdge(graph, 6, 7);
    addEdge(graph, 7, 5);
    addEdge(graph, 3, 7);
    addEdge(graph, 4,6);
    addEdge(graph, 1,5);
    addEdge(graph, 2,4);
    addEdge(graph, 0,3);
    addEdge(graph, 0,4);
    addEdge(graph, 0,7);
    */
    
    
    
    
    
    
 
    return 0;
}
