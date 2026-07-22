#include "dijkstras.h"
#include <stdlib.h>

typedef struct
{
    GraphScalarType minCost;
    GraphNode *prevNode;
} ShortestPathState;

// Heap
typedef struct
{
    GraphScalarType weight;
    int nodeIndex;
    int fromIndex;
} HeapElem;

typedef struct
{
    HeapElem *arr;
    int count;
    int capacity;
} Heap;

Heap *create_heap(int capacity);
void heap_insert(Heap *h, HeapElem key);
void heapify_bottom_top(Heap *h, int index);
void heapify_top_bottom(Heap *h, int parent_node);
HeapElem heap_pop_min(Heap *h);

#define REAL_HEAP 1

Heap *create_heap(int capacity)
{
    Heap *h = (Heap *)malloc(sizeof(Heap));
    if (h == NULL)
    {
        return NULL;
    }

    h->count = 0;
    h->capacity = capacity;

    h->arr = (HeapElem *)malloc(capacity * sizeof(HeapElem));
    if (h->arr == NULL)
    {
        return NULL;
    }

    return h;
}

#ifdef REAL_HEAP

__attribute__((always_inline)) void heap_insert(Heap *h, HeapElem elem)
{
    h->arr[h->count] = elem;
    heapify_bottom_top(h, h->count);
    h->count++;
}

__attribute__((always_inline)) void heapify_bottom_top(Heap *h, int index)
{
    HeapElem temp;
    int parent_node = (index - 1) / 2;
    int tempIndex;

    while (h->arr[parent_node].weight > h->arr[index].weight)
    {
        // swap and recursive call
        temp = h->arr[parent_node];
        h->arr[parent_node] = h->arr[index];
        h->arr[index] = temp;
        tempIndex = index;
        parent_node = (index - 1) / 2;
        index = tempIndex;
    }
}

__attribute__((always_inline)) void heapify_top_bottom(Heap *h,
                                                       int parent_node)
{
    bool cont;
    do
    {
        int left = parent_node * 2 + 1;
        int right = parent_node * 2 + 2;
        HeapElem temp;
        int min;

        if (left >= h->count || left < 0)
            left = -1;
        if (right >= h->count || right < 0)
            right = -1;

        if (left != -1 && h->arr[left].weight < h->arr[parent_node].weight)
            min = left;
        else
            min = parent_node;

        if (right != -1 && h->arr[right].weight < h->arr[min].weight)
            min = right;

        cont = false;
        if (min != parent_node)
        {
            temp = h->arr[min];
            h->arr[min] = h->arr[parent_node];
            h->arr[parent_node] = temp;
            parent_node = min;
            cont = true;
        }
    } while (cont);
}

__attribute__((always_inline)) HeapElem heap_pop_min(Heap *h)
{
    HeapElem pop;
    // replace first node by last and delete last
    pop = h->arr[0];
    h->arr[0] = h->arr[h->count - 1];
    h->count--;
    heapify_top_bottom(h, 0);
    return pop;
}

#else

void heap_insert(Heap *h, HeapElem key)
{
    h->arr[h->count] = key;
    h->count++;
}

HeapElem heap_pop_min(Heap *h)
{
    int minIdx = -1;
    int minVal = 0x7FFFFFFF;

    for (int i = 0; i < h->count; i++)
    {
        if (h->arr[i].weight < minVal)
        {
            minVal = h->arr[i].weight;
            minIdx = i;
        }
    }

    HeapElem ret = h->arr[minIdx];

    for (int i = minIdx + 1; i < h->count; i++)
    {
        h->arr[i - 1] = h->arr[i];
    }

    h->count--;

    return ret;
}

#endif

__efficient__ void shortest_path_internal(Heap *h, ShortestPathState *state,
                                          Graph *g, int to)
{
    while (h->count)
    {
        HeapElem curr = heap_pop_min(h);

        if (state[curr.nodeIndex].minCost != INF_VALUE)
        {
            continue;
        }

        state[curr.nodeIndex].minCost = curr.weight;
        state[curr.nodeIndex].prevNode = &g->nodes[curr.fromIndex];

        if (curr.nodeIndex == to)
        {
            break;
        }

        GraphNode node = g->nodes[curr.nodeIndex];
        for (int i = 0; i < node.adjLen; i++)
        {
            HeapElem e;
            e.fromIndex = curr.nodeIndex;
            e.nodeIndex = node.adjList[i].destIdx;
            e.weight = node.adjList[i].cost + curr.weight;
            heap_insert(h, e);
        }
    }
}

GraphEdge *shortest_path(Graph *g, GraphNode *from, GraphNode *to,
                         int *pathLength)
{
    // Creating and initializing the state
    ShortestPathState *state =
        (ShortestPathState *)malloc(sizeof(ShortestPathState) * g->nodeLen);
    if (!state)
    {
        return NULL;
    }

    int edgeCount = 0;

    for (int i = 0; i < g->nodeLen; i++)
    {
        state[i].minCost = INF_VALUE;
        edgeCount += g->nodes[i].adjLen;
    }

    Heap *minHeap = create_heap(edgeCount + 1);
    HeapElem root;
    root.weight = 0;
    root.nodeIndex = from->index;
    root.fromIndex = -1;
    heap_insert(minHeap, root);

    shortest_path_internal(minHeap, state, g, to->index);

    // freeing heap
    free(minHeap->arr);
    free(minHeap);

    if (state[to->index].minCost == INF_VALUE)
    {
        free(state);
        return NULL;
    }

    // finding the path
    // first, obtaining the path length
    int currIdx = to->index;
    int length = 0;
    while (currIdx != from->index)
    {
        length++;
        currIdx = state[currIdx].prevNode->index;
    }

    // allocating the path and storing the edges
    GraphEdge *path = (GraphEdge *)malloc(sizeof(GraphEdge) * length);
    if (!path)
    {
        free(state);
        return NULL;
    }

    currIdx = to->index;
    int counter = length - 1;
    while (currIdx != from->index)
    {
        GraphEdge e;
        e.destIdx = currIdx;
        e.cost = state[currIdx].minCost -
                 state[state[currIdx].prevNode->index].minCost;
        path[counter] = e;
        currIdx = state[currIdx].prevNode->index;
        counter--;
    }

    *pathLength = length;

    free(state);

    return path;
}
