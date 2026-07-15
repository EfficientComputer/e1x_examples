#include <stdio.h>
#include "dijkstras.h"
#include <stdlib.h>

#define NUM_ITERATIONS 1

#define MAP_WIDTH 20
#define MAP_HEIGHT 9

unsigned char map[] =
    "\
S****+****+*****+ + \
     *    *     * * \
+****+*+ ++**+  +*+ \
*            *      \
+***+**+  D  +****+ \
    *     *       * \
    +*****+*******+ \
    *             + \
    +               \
";

typedef struct
{
    int x;
    int y;
} VertexData;

void load_graph_helper(Graph *g, unsigned char *map, int x, int y, int src,
                       int count, int *start, int *dest)
{
    if (x < 0 || x >= MAP_WIDTH || y < 0 || y >= MAP_HEIGHT)
    {
        return;
    }

    unsigned char v = map[MAP_WIDTH * y + x];

    if (v == ' ')
    {
        return;
    }
    else if (v == 'S' || v == 'D' || v == '+')
    {
        // new vertex
        VertexData *data = (VertexData *)malloc(sizeof(VertexData));
        data->x = x;
        data->y = y;

        int nodeIdx = graph_add_node(g, data);
        map[MAP_WIDTH * y + x] = 128 + nodeIdx;

        if (src >= 0)
        {
            graph_node_add_edge(&g->nodes[src], nodeIdx, count);
            graph_node_add_edge(&g->nodes[nodeIdx], src, count);
        }

        count = 0;
        src = nodeIdx;

        if (v == 'S')
        {
            *start = nodeIdx;
        }
        else if (v == 'D')
        {
            *dest = nodeIdx;
        }
    }
    else if (v >= 128)
    {
        // existing node - must make sure that it doesn't already have an edge
        // between src and existing
        GraphNode *existing = &g->nodes[v - 128];
        bool alreadyConnected = false;
        for (int i = 0; i < existing->adjLen; i++)
        {
            if (existing->adjList[i].destIdx == src)
            {
                alreadyConnected = true;
                break;
            }
        }

        if (!alreadyConnected && existing->index != src)
        {
            int existingIndex = existing->index;
            graph_node_add_edge(existing, src, count);
            graph_node_add_edge(&g->nodes[src], existingIndex, count);
        }

        return;
    }
    else
    {
        map[MAP_WIDTH * y + x] = ' ';
    }

    load_graph_helper(g, map, x, y + 1, src, count + 1, start, dest);
    load_graph_helper(g, map, x, y - 1, src, count + 1, start, dest);
    load_graph_helper(g, map, x + 1, y, src, count + 1, start, dest);
    load_graph_helper(g, map, x - 1, y, src, count + 1, start, dest);
}

void load_graph(Graph *g, unsigned char *map, int *start, int *dest)
{
    load_graph_helper(g, map, 0, 0, -1, 0, start, dest);
}

GraphEdge *shortest_path(Graph *g, GraphNode *from, GraphNode *to,
                         int *pathLength);

void warmup()
{
    Graph *g = new_graph();
    int startIdx = graph_add_node(g, NULL);
    int dstIdx = graph_add_node(g, NULL);
    graph_node_add_edge(&g->nodes[startIdx], dstIdx, 1);

    int pathLength;
    GraphEdge *path =
        shortest_path(g, &g->nodes[startIdx], &g->nodes[dstIdx], &pathLength);

    free(path);

    // deleting the graph
    int gLength = g->nodeLen;
    for (int i = gLength - 1; i >= 0; i--)
    {
        void *data = graph_delete_node(g, i);
    }

    free(g->nodes);
    free(g);
}

int main()
{
    // warmup();

    Graph *g = new_graph();
    int start;
    int dest;
    load_graph(g, map, &start, &dest);

    int pathLength;
    GraphEdge *path = NULL;
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
    {
        if (path != NULL)
            free(path);
        path = shortest_path(g, &g->nodes[start], &g->nodes[dest], &pathLength);
    }

    // printing the path
    VertexData *startData = (VertexData *)g->nodes[start].data;
    printf("START: (%d, %d)\n", startData->x, startData->y);
    GraphScalarType distance = 0;
    for (int i = 0; i < pathLength; i++)
    {
        VertexData *data = (VertexData *)g->nodes[path[i].destIdx].data;
        printf("    -> (%d, %d)\n", data->x, data->y);
        distance += path[i].cost;
    }
    printf("PATH DISTANCE: %u\n", distance);

    free(path);

    // deleting the graph
    int gLength = g->nodeLen;
    for (int i = gLength - 1; i >= 0; i--)
    {
        void *data = graph_delete_node(g, i);
        free(data);
    }

    free(g->nodes);
    free(g);

    if (pathLength != 8)
    {
        printf(
            "[dijkstras] FAIL - Path segment count of %d did not match "
            "expected "
            "(8)\n",
            pathLength);
        return 1;
    }
    else if (distance != 28)
    {
        printf(
            "[dijkstras] FAIL - Path distance of %d did not match expected "
            "(28)\n",
            distance);
        return 1;
    }

    printf("[dijkstras] PASS\n");

    return 0;
}
