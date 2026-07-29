#define BFS_NUM_VTX 1000
#define BFS_NUM_EDGE 10000

#ifdef EFF_BLD_HAND_OPTIMIZED
void bfs(unsigned startv, const unsigned *const restrict oa,
         const unsigned *const restrict na, unsigned *frontier);
#else
void bfs(unsigned start_vertex, const unsigned *const vertices,
         const unsigned *const edge_destinations, unsigned *frontier);

#endif
