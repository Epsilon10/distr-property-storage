#include "galois/Galois.h"
#include "galois/graphs/LS_LC_CSR_64_Graph.h"

#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <atomic>
#include <vector>
#include <thread>
#include <mutex>
#include <initializer_list>
#include <cassert>
#include <chrono>
#include <memory>

#include "galois/Bag.h"
#include "galois/Reduction.h"

#include "DistBench/Output.h"
#include "DistBench/Start.h"
#include "galois/DistGalois.h"
#include "galois/gstl.h"
#include "galois/DReducible.h"
#include "galois/DTerminationDetector.h"
#include "galois/runtime/Tracer.h"

#include "queue.h"

#include "buckets.h"

BucketContainer* buckets;

#include "ds_sync_struct.h"

constexpr uint64_t NUM_SLOTS = 1 << 8;
constexpr uint64_t NUM_NODES = 1 << 20;

enum RequestType {
    LIGHT, HEAVY, BOTH
};

const uint64_t infinity = std::numeric_limits<uint64_t>::max();

using Graph = galois::graphs::DistGraph<NodeData, uint64_t>;

std::unique_ptr<galois::graphs::GluonSubstrate<Graph>> syncSubstrate;
galois::DynamicBitSet bitset_dist_current;

struct InitializeGraph {
    Graph* graph;
    std::uint64_t local_src;

    InitializeGraph(Graph* graph, std::uint64_t local_src) 
        : graph(graph), local_src(local_src) {
    }

    void operator()(std::uint64_t src) const {
        NodeData& sdata = graph->getData(src);
        sdata.dist = src == local_src ? 0 : infinity;
        //std::cout << "sdata " << sdata << std::endl;
    }
};

GALOIS_SYNC_STRUCTURE_BITSET(dist_current);

void delta_step(Graph& graph, std::uint64_t src, BucketContainer* buckets);

//////////
// MAIN //
//////////
constexpr static const char* const name =
    "SSSP - Distributed Heterogeneous with "
    "worklist.";
constexpr static const char* const desc = "BFS on Distributed Galois.";
constexpr static const char* const url  = nullptr;

int main(int argc, char** argv) {
    // TODO: add serialzation code
    galois::DistMemSys G;
    DistBenchStart(argc, argv, name, desc, url);

    uint64_t src_node = 0;
    std::srand(std::time(NULL));

    const auto& net = galois::runtime::getSystemNetworkInterface();

    std::unique_ptr<Graph> hg;
    std::tie(hg, syncSubstrate) = distGraphInitialization<NodeData, uint64_t>();

    buckets = new BucketContainer(NUM_SLOTS);

   // Graph graph = load_csr_64_graph(); // does parallel ingestion
    bitset_dist_current.resize(hg->size());

   const auto& allNodes = hg->allNodesRange();
   
    galois::do_all(
        galois::iterate(allNodes.begin(), allNodes.end()),
        [&](std::uint64_t i) {
            if (i == 0) {
                hg->getData(i) = 0;
            } else {
                hg->getData(i) = infinity;
            }
        });

    auto start = std::chrono::high_resolution_clock::now();
    delta_step(*hg, 0, buckets);
    auto end = std::chrono::high_resolution_clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "TIME: " << time_ms.count() << std::endl;

    for (int i = 0; i < NUM_NODES; i++) {
        uint64_t dist = hg->getData(i).dist;
        std::cout << "VERTEX: " << i << "   FINAL DIST " << hg->getData(i).dist << std::endl;
    }
}

bool relax_edges_for_vertex(
    Graph& graph, BucketContainer* buckets, std::uint64_t u, RequestType req_type, 
galois::GReduceMin<std::uint64_t>& distributed_min, std::uint64_t cur_bucket, Queue<BucketContainer::Vertex>& same_bucket) {
    std::atomic<bool> any_updated = false;
    auto degree = graph.getDegree(u);
    auto edge_itr_start = graph.edge_begin(u);
    auto const u_dist = graph.getData(u, galois::MethodFlag::READ);
    auto edge_itr_end = graph.edge_end(u); // don't need this
    std::mutex same_bucket_lock;
    
    galois::do_all(galois::iterate((uint64_t)0, (uint64_t)degree), 
        [&](std::uint64_t i) {
            // iterate through each neighbor v
            auto v = graph.getEdgeDst(edge_itr_start + i);

            NodeData& node_data = graph.getData(v, galois::MethodFlag::WRITE);
            std::uint64_t cur_dist = node_data.dist;
            // need write perms because we may relax v's current distance
            auto const& w = graph.getEdgeData(u, v); // dist between u and v

            std::uint64_t new_dist = w + u_dist;

            if (LIGHT || BOTH) {
                if (w <= DELTA) {
                    
                    // relax the edge
                    if (new_dist < cur_dist) {
                        if (!any_updated.load(std::memory_order_relaxed))
                            any_updated.store(true, std::memory_order_relaxed);
                        

                        // always update distance, even if not local
                        // we set the corresponding bit and it is communicated 
                        cur_dist = new_dist;
                        uint64_t new_bucket = cur_dist / DELTA;
                        bitset_dist_current.set(v);

                        // only update buckets if local
                       if (graph.isOwned(v)) {
                            if (new_bucket != cur_bucket)
                                distributed_min.update(new_bucket);
                            if (new_bucket == cur_bucket && req_type == LIGHT) {
                                same_bucket_lock.lock();
                                same_bucket.add(new BucketContainer::Vertex(cur_bucket));
                                same_bucket_lock.unlock();
                            } else {
                            //   std::cout << "GETTING BUCKET " << new_bucket << std::endl;
                                buckets->push(new_bucket, v);
                            }
                       }
                    }
                }
            }

            if (HEAVY || BOTH) {
                if (w > DELTA) {
                    // relax the edge
                    if (new_dist < cur_dist) {
                        if (!any_updated.load(std::memory_order_relaxed))
                            any_updated.store(true, std::memory_order_relaxed);

                        cur_dist = new_dist;

                        bitset_dist_current.set(v);
                        uint64_t new_bucket = cur_dist / DELTA;

                        if (graph.isOwned(v)) {
                            if (new_bucket != cur_bucket)
                                distributed_min.update(new_bucket);

                            if (new_bucket == cur_bucket && req_type == LIGHT) {
                                same_bucket_lock.lock();
                                same_bucket.add(new BucketContainer::Vertex(cur_bucket));
                                same_bucket_lock.unlock();
                            } else {
                                buckets->push(new_bucket, v);
                            }
                        }
                    }
                }
            }
            
    },
    galois::loopname("relax_edges"), galois::steal()
    );

    syncSubstrate->sync<writeDestination, readSource, Reduce_min_bucket,
                        Bitset_dist_current, false>("poop");
    
    //std::cout << "made it out" << std::endl;
    return any_updated;
 }

 void delta_step(Graph& graph, std::uint64_t src, BucketContainer* buckets) {
    galois::GReduceMin<std::uint64_t> distributed_min; // starts with uint64_t max
    std::uint64_t min_bucket = 0; // find intl
    Queue<BucketContainer::Vertex> dummy_bag;

    relax_edges_for_vertex(graph,buckets,0, BOTH, distributed_min, min_bucket, dummy_bag);
    
    if (buckets->is_empty(min_bucket)) {
        min_bucket = distributed_min.reduce();
    }

    distributed_min.reset();
    std::cout << "LESS GET STARTED\n";
    // only relax the minimum bucket

    while (true)  { // fix cond
        Bucket* bucket = buckets->get(min_bucket); // retrieves min insert bag, ref
        std::unordered_set<std::uint64_t> r;
        
        bool any_updated_light;
        
        while (!bucket.empty()) {
            uint64_t cnt = 0;

            Queue<BucketContainer::Vertex> cur_bucket_replace;

            bucket->vertices.all([&graph, buckets, &distributed_min, min_bucket, cur_bucket_replace, &r](BucketContainer::Vertex* vtx) { 
                uint64_t vertex = vtx->val;
                any_updated_light = relax_edges_for_vertex(graph,buckets, vertex,LIGHT, distributed_min, min_bucket, cur_bucket_replace);
                r.insert(vertex);
            }); // 
            bucket->vertices.remove_all(); // clears old contents
            bucket->vertices = std::move(cur_bucket_replace);
        }

        // now the bucket is empty, we can free the assoc heap memory
        buckets->remove_bucket(min_bucket);

       // std::cout << "ROUND 2\n";
        bool any_updated_heavy;
        for (std::uint64_t vertex : r) {
            any_updated_heavy = relax_edges_for_vertex(graph,buckets, vertex,HEAVY, distributed_min, min_bucket, dummy_bag); // will; now refill bucket 
        }
        if (!any_updated_heavy && !any_updated_light) {
            std::cout << "FINAL STEP\n";
            uint64_t test = ~0U;
            for (uint64_t i = 0; i < NUM_BUCKETS; i++) {
                if (!buckets->is_empty(i)) {
                    test = std::min(test, buckets->get_min(i));
                }
            }
            if (test == ~0U) break;
            std::cout << "TEST VAL " << test << std::endl;
            min_bucket = test;
        } else {
            min_bucket = distributed_min.reduce();
        // std::cout << "UPDATDE MIN BUCKET: " << min_bucket << std::endl;
        }
        distributed_min.reset();
        std::cout << "UPDATED MIN BUCKET " << min_bucket << std::endl;
       // std::cout << "VERTEX 662: " << graph.getData(662) << std::endl;
    }
 }  
