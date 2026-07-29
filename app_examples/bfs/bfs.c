#include <stdint.h>

#include "bfs.h"

#ifdef EFF_BLD_HAND_OPTIMIZED
extern int __effcc_nop(int);

static inline int visit_vtx(int v, const unsigned *const restrict oa,
                            const unsigned *const restrict na,
                            unsigned *const restrict frontier, int iter,
                            int next)
{
    int v_buf = __effcc_nop(v);
    int tv = __atomic_load_n(frontier + v_buf, __ATOMIC_ACQUIRE);
    __atomic_store_n(frontier + v_buf, tv, __ATOMIC_RELEASE);
    unsigned new_tv = __effcc_nop(tv ^ v_buf) ^ v_buf;

    unsigned start_offset = oa[v];
    unsigned end_offset = oa[v + 1];
    unsigned ch = 0;

    if (new_tv == iter)
    {
        for (int i = start_offset; i < end_offset; i++)
        {
            unsigned nei = na[i];

            unsigned nei_buf = __effcc_nop(nei);
            unsigned *frontier_nei = frontier + nei_buf;
            unsigned tn = __atomic_load_n(frontier + nei_buf, __ATOMIC_ACQUIRE);
            unsigned newval = (tn == BFS_NUM_VTX) ? next : tn;
            __atomic_store_n(frontier + nei_buf, newval, __ATOMIC_RELEASE);
            unsigned new_tn = __effcc_nop(tn ^ nei_buf) ^ nei_buf;
            ch |= (new_tn == BFS_NUM_VTX);
        }
    }
    return ch;
}

__efficient__ void bfs(unsigned startv, const unsigned *const restrict oa,
                       const unsigned *const restrict na,
                       unsigned *const restrict frontier)
{
    /*This function runs a BFS traversal over a graph stored in CSR format.
      OA is the Offsets Array for the CSR and NA is the Neighbors Array for
      the CSR.  OA must contain #Vertices + 1 entries and NA contains one
      entry per vertex.

      The algorithm is written to favor parallelism over work efficiency,
      opting for a parallel scan of all entries in the frontier, to find
      active elements, rather than a more dynamic data structure that addsC
      potential impediments to parallelism.

      The output in the frontier array is the discovery level of each
      vertex.

    */

    unsigned not_done = 1;
    unsigned iter = 0;
    frontier[startv] = 0;

    while (not_done)
    {
        not_done = 0;
        unsigned next = iter + 1;

        /*
          These each execute in parallel over sub-ranges of the vertices
          For now, this is written to be 6x parallel, which seems to fill
          the e1x fabric in nicely.
        */

        int ch = 0;
        for (int v = 0; v < BFS_NUM_VTX / 6; v++)
        {
            ch |= visit_vtx(v, oa, na, frontier, iter, next);
        }

        int ch2 = 0;
        for (int v = BFS_NUM_VTX / 6; v < BFS_NUM_VTX * 2 / 6; v++)
        {
            ch2 |= visit_vtx(v, oa, na, frontier, iter, next);
        }

        int ch3 = 0;
        for (int v = BFS_NUM_VTX * 2 / 6; v < BFS_NUM_VTX * 3 / 6; v++)
        {
            ch3 |= visit_vtx(v, oa, na, frontier, iter, next);
        }

        int ch4 = 0;
        for (int v = BFS_NUM_VTX * 3 / 6; v < BFS_NUM_VTX * 4 / 6; v++)
        {
            ch4 |= visit_vtx(v, oa, na, frontier, iter, next);
        }

        int ch5 = 0;
        for (int v = BFS_NUM_VTX * 4 / 6; v < BFS_NUM_VTX * 5 / 6; v++)
        {
            ch5 |= visit_vtx(v, oa, na, frontier, iter, next);
        }
        int ch6 = 0;
        for (int v = BFS_NUM_VTX * 5 / 6; v < BFS_NUM_VTX; v++)
        {
            ch6 |= visit_vtx(v, oa, na, frontier, iter, next);
        }

        not_done |= (ch | ch2 | ch3 | ch4 | ch5 | ch6);

        iter++;
    }
}
#else
__efficient__ void bfs(
    unsigned start_vertex,                   // the index of the starting vertex
    const unsigned *const vertices,          // the array of offsets into the
                                             // edge_destinations array for each vertex
    const unsigned *const edge_destinations, // the array of destination vertex
                                             // indicies for each edge
    unsigned *frontier)                      // the array of distances from the start vertex
{
    const unsigned vertex_count = BFS_NUM_VTX;
    const unsigned edge_count = BFS_NUM_EDGE;

    static uint32_t visited[BFS_NUM_VTX]; // the array of visited flags for each vertex
    static int32_t queue[BFS_NUM_VTX];    // vertex index queue

    for (int i = 0; i < BFS_NUM_VTX; i++)
    {
        visited[i] = 0;
        queue[i] = 0;
    }

    queue[0] = start_vertex;
    frontier[start_vertex] = 0;
    int queue_pointer = 0;
    int queue_length = 1;
    visited[start_vertex] = 1;

    while (queue_length > 0)
    {
        int32_t vertex_idx = queue[queue_pointer];
        queue_length--;

        // incrementing the queue pointer and wrapping around if necessary
        queue_pointer++;
        if (queue_pointer >= vertex_count)
        {
            queue_pointer -= vertex_count;
        }

        uint32_t start_offset = vertices[vertex_idx];
        uint32_t end_offset = edge_count;
        if (vertex_idx != vertex_count - 1)
        {
            end_offset = vertices[vertex_idx + 1];
        }

        for (int i = start_offset; i < end_offset; i++)
        {
            uint32_t neighbor_idx = edge_destinations[i];

            if (visited[neighbor_idx] == 0)
            {
                int32_t next_queue_pointer = queue_pointer + queue_length;
                if (next_queue_pointer >= vertex_count)
                {
                    next_queue_pointer -= vertex_count;
                }

                queue[next_queue_pointer] = neighbor_idx;
                visited[neighbor_idx] = 1;
                queue_length++;
                frontier[neighbor_idx] = frontier[vertex_idx] + 1;
            }
        }
    }
}
#endif
