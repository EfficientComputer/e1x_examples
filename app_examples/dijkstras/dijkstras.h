#include <stdbool.h>
#include <stdint.h>

typedef uint32_t GraphScalarType;
#define INF_VALUE 0xFFFFFFFF

struct GraphNode;

typedef struct
{
    int destIdx;
    GraphScalarType cost;
} GraphEdge;

typedef struct GraphNode
{
    GraphEdge *adjList;
    void *data;
    int adjLen;
    int adjCap;
    int index;
} GraphNode;

typedef struct
{
    GraphNode *nodes;
    int nodeLen;
    int nodeCap;
} Graph;

Graph *new_graph();

int graph_add_node(Graph *g, void *data);
void *graph_delete_node(Graph *g, int index);
bool graph_node_add_edge(GraphNode *n, int destIdx, GraphScalarType cost);
bool graph_node_remove_edge(GraphNode *n, int destIdx);
