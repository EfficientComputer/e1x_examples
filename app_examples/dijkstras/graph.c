#include "dijkstras.h"
#include <stdlib.h>

#define ADJ_LIST_DEFAULT_CAP 4
#define NODES_DEFAULT_CAP 16

Graph *new_graph()
{
    Graph *ret = (Graph *)malloc(sizeof(Graph));
    if (!ret)
    {
        return NULL;
    }

    ret->nodeCap = NODES_DEFAULT_CAP;
    ret->nodeLen = 0;
    ret->nodes = (GraphNode *)malloc(sizeof(GraphNode) * NODES_DEFAULT_CAP);

    return ret;
}

bool graph_resize_node_arr(Graph *g, int newCap)
{
    GraphNode *newArr =
        (GraphNode *)realloc(g->nodes, sizeof(GraphNode) * newCap);

    if (!newArr)
    {
        return false;
    }

    // updating graph state
    g->nodeCap = newCap;
    g->nodes = newArr;

    return true;
}

int graph_add_node(Graph *g, void *data)
{
    if (g->nodeCap == g->nodeLen)
    {
        if (!graph_resize_node_arr(g, 2 * g->nodeCap))
        {
            return -1;
        }
    }

    GraphNode n;
    n.data = data;
    n.index = g->nodeLen;
    n.adjCap = ADJ_LIST_DEFAULT_CAP;
    n.adjLen = 0;
    n.adjList = (GraphEdge *)malloc(sizeof(GraphEdge) * ADJ_LIST_DEFAULT_CAP);

    g->nodes[g->nodeLen] = n;
    g->nodeLen++;
    return g->nodeLen - 1;
}

void *graph_delete_node(Graph *g, int index)
{
    free(g->nodes[index].adjList);

    void *data = g->nodes[index].data;

    // shifting all the existing nodes down
    for (int i = 1; i < g->nodeLen; i++)
    {
        if (i > index)
        {
            g->nodes[i - 1] = g->nodes[i];
            g->nodes[i - 1].index = i - 1;
        }

        // updating edges
        for (int j = 0; j < g->nodes[i - 1].adjLen; j++)
        {
            int edgeDestIdx = g->nodes[i - 1].adjList[j].destIdx;
            if (edgeDestIdx == index)
            {
                graph_node_remove_edge(&g->nodes[i - 1], index);
            }
            else if (edgeDestIdx > index)
            {
                g->nodes[i - 1].adjList[j].destIdx--;
            }
        }
    }

    g->nodeLen--;

    // if length is less than a third the capacity, shrink the array
    if (g->nodeCap > 3 * g->nodeLen && g->nodeCap > NODES_DEFAULT_CAP)
    {
        graph_resize_node_arr(g, g->nodeCap / 2);
    }

    return data;
}

bool graph_node_resize_adj_arr(GraphNode *n, int newCap)
{
    // grows capacity by 2x each time
    GraphEdge *newArr =
        (GraphEdge *)realloc(n->adjList, sizeof(GraphEdge) * newCap);
    if (!newArr)
    {
        return false;
    }

    // updating graph state
    n->adjCap = newCap;
    n->adjList = newArr;
    return true;
}

bool graph_node_add_edge(GraphNode *n, int destIdx, GraphScalarType cost)
{
    if (n->adjCap == n->adjLen)
    {
        if (!graph_node_resize_adj_arr(n, 2 * n->adjCap))
        {
            return false;
        }
    }

    GraphEdge e;
    e.destIdx = destIdx;
    e.cost = cost;

    n->adjList[n->adjLen] = e;
    n->adjLen++;

    return true;
}

bool graph_node_remove_edge(GraphNode *n, int destIdx)
{
    bool shifting = false;
    for (int i = 0; i < n->adjLen; i++)
    {
        if (!shifting && n->adjList[i].destIdx == destIdx)
        {
            shifting = true;
        }
        else if (shifting)
        {
            n->adjList[i - 1] = n->adjList[i];
        }
    }

    if (shifting)
    {
        n->adjLen--;
    }

    // if length is less than a third the capacity, shrink the array
    if (n->adjCap > 3 * n->adjLen)
    {
        if (!graph_node_resize_adj_arr(n, n->adjCap / 2))
        {
            return false;
        }
    }

    return true;
}
