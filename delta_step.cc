#include "galois/Galois.h"
#include "galois/graphs/LS_LC_CSR_64_Graph.h"
#include "loader.hpp"

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

#include "galois/Bag.h"
#include "galois/Reduction.h"

std::uint64_t constexpr DELTA = 80;

std::uint64_t constexpr NUM_BUCKETS = 1 << 8;

enum RequestType {
    LIGHT, HEAVY, BOTH
};

const uint64_t infinity = std::numeric_limits<uint64_t>::max();

using Graph = galois::graphs::Ben_CSR_64<std::uint64_t, EdgeData_N, true>;

template<typename T>
struct ConcurrentSparseTable {
    private:
    struct Elem {
        std::uint64_t tag;
        T val;
        Elem* next = nullptr;
        Elem() = default;
        ~Elem() = default;
    };

    std::uint64_t capacity; 

    std::vector<Elem*> raw_array;
    std::vector<std::mutex> locks;
    std::vector<std::uint64_t> min_bucket_arr;

    public:
    ConcurrentSparseTable(std::size_t capacity) : capacity(capacity), raw_array(capacity, nullptr), locks(capacity), min_bucket_arr(capacity, ~0U)  {
        
    }
    
    public:
    bool is_empty(std::uint64_t i) {
        return raw_array[i] == nullptr;
    }
    
    // WARNING ONLY CALL WITHIN LOCKED REGION
    void update_min(uint64_t idx) {
        Elem* cur = raw_array[idx];
        uint64_t min_tag = ~0U;

        while (cur != nullptr) {
            min_tag = std::min(cur->tag, min_tag);
            cur = cur->next;
        }

        min_bucket_arr[idx] = min_tag;
    }

    template<typename Func>
    void push(std::uint64_t i, std::uint64_t v, Func f) {
        std::uint64_t idx = i & (capacity - 1);
     //   std::cout << "Acquring lock " << idx << " for indice " << i << std::endl;
        std::lock_guard g{locks[idx]};
     //   std::cout << "hi\n";
        
        Elem* node = raw_array[idx];

        if (node == nullptr) {
            Elem* x = new Elem();
            x->tag = i;
            raw_array[idx] = x;
            min_bucket_arr[idx] = i;
            f(x->val);
            return;
        }

        Elem* last = nullptr;
        while (node != nullptr) {
            if (node->tag == i) {
                f(node->val);
                return;
            }
            last = node;
            node = node->next;
        }

        // tag not there
        Elem* to_add = new Elem();
        to_add->tag = i;
        last->next = to_add;

        min_bucket_arr[idx] = std::min(min_bucket_arr[idx], i);
        f(to_add->val);
    }

    T& operator[](std::uint64_t i) {
        std::uint64_t idx = i & (capacity - 1);
     //   std::cout << "Acquring lock " << idx << " for indice " << i << std::endl;
     //   std::cout << "hi\n";
        
        Elem* node = raw_array[idx];

        if (node == nullptr) {
            Elem* x = new Elem();
            x->tag = i;
            raw_array[idx] = x;
            return x->val;
        }

        Elem* last = nullptr;
        while (node != nullptr) {
            if (node->tag == i)
                return node->val;
            last = node;
            node = node->next;
        }

        // tag not there
        Elem* to_add = new Elem();
        to_add->tag = i;
        last->next = to_add;

        return to_add->val;
    }

    void erase(std::uint64_t i) {
        std::uint64_t idx = i & (capacity - 1);
        std::lock_guard g{locks[idx]};
       // std::cout << "erasing " << i << std::endl;

        Elem* node = raw_array[idx];
        assert(node != nullptr);

        if (node->next == nullptr) {
            assert(node->tag == i);
            delete node;
            raw_array[idx] = nullptr;            
            return;
        }

        if (node->tag == i) {
            raw_array[idx] = node->next;
            delete node;
            update_min(idx);
            return;
        }

        Elem* last = nullptr;
        while (node != nullptr) {
            if (node->tag == i) {
                last->next = node->next;
                delete node;
                break;
            }
            last = node;
            node = node->next;
        }
        update_min(idx);
    }

    uint64_t get_min(uint64_t i) {
        return min_bucket_arr[i];
    }
};

using BucketContainer = ConcurrentSparseTable<galois::InsertBag<std::uint64_t>>;

struct InitializeGraph {
    Graph* graph;
    std::uint64_t local_src;

    InitializeGraph(Graph* graph, std::uint64_t local_src) 
        : graph(graph), local_src(local_src) {
    }

    void operator()(std::uint64_t src) const {
        std::uint64_t& sdata = graph->getData(src);
        sdata = src == local_src ? 0 : infinity;
        //std::cout << "sdata " << sdata << std::endl;
    }
};

void delta_step(Graph& graph, std::uint64_t src, BucketContainer& buckets);

int main() {
    // TODO: add serialzation code
    galois::SharedMemSys G;
    uint64_t src_node = 0;
    std::srand(std::time(NULL));

    BucketContainer buckets{256};

    Graph graph = load_csr_64_graph(); // does parallel ingestion
   // std::cout << "hi\n";
    galois::do_all(
        galois::iterate((uint64_t)0, NUM_NODES),
        [&](std::uint64_t i) {
            if (i == 0) {
                graph.getData(i) = 0;
            } else {
                graph.getData(i) = infinity;
               // std::cout << graph.getData(i) << std::endl;
            }
        });

  //  std::cout << "hi\n";
    auto start = std::chrono::high_resolution_clock::now();
    delta_step(graph, 0, buckets);
    auto end = std::chrono::high_resolution_clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "TIME: " << time_ms.count() << std::endl;

    for (int i = 0; i < NUM_NODES; i++) {
        uint64_t dist = graph.getData(i);
        std::cout << "VERTEX: " << i << "   FINAL DIST " << graph.getData(i) << std::endl;
    }
}

bool relax_edges_for_vertex(
    Graph& graph, BucketContainer& buckets, std::uint64_t u, RequestType req_type, 
    galois::GReduceMin<std::uint64_t>& distributed_min, std::uint64_t cur_bucket, galois::InsertBag<std::uint64_t>& sbb) {
    std::atomic<bool> any_updated = false;
    auto degree = graph.getDegree(u);
    auto edge_itr_start = graph.edge_begin(u);
    auto const u_dist = graph.getData(u, galois::MethodFlag::READ);
    auto edge_itr_end = graph.edge_end(u); // don't need this
    
    galois::do_all(galois::iterate((uint64_t)0, (uint64_t)degree), 
        [&](std::uint64_t i) {
            // iterate through each neighbor v
            auto v = graph.getEdgeDst(edge_itr_start + i);
           // std::cout << "u: " << u << std::endl;
           // std::cout << "v: " << v << std::endl;
            std::uint64_t& cur_dist = graph.getData(v, galois::MethodFlag::WRITE);
            // need write perms because we may relax v's current distance
            auto const& w = graph.getEdgeData(u, v); // dist between u and v
          //  std::cout << "W: " << w.weight << std::endl;
            std::uint64_t new_dist = w + u_dist;
           // std::cout << "new dist " << new_dist << std::endl;
            if (LIGHT || BOTH) {
                if (w <= DELTA) {
                    
                    // relax the edge
                    if (new_dist < cur_dist) {
                        if (!any_updated.load(std::memory_order_relaxed))
                            any_updated.store(true, std::memory_order_relaxed);
                       // std::cout << "cur dist: " << cur_dist << "    new d " << new_dist << std::endl;
                        cur_dist = new_dist;
                        uint64_t new_bucket = cur_dist / DELTA;
                            //                    std::cout << "new d " << new_dist << std::endl;

                       // std::cout << "new bucket: " << new_bucket << std::endl;
                       if (new_bucket != cur_bucket)
                            distributed_min.update(new_bucket);
                        if (new_bucket == cur_bucket && req_type == LIGHT) {
                            sbb.push(v);
                        } else {
                         //   std::cout << "GETTING BUCKET " << new_bucket << std::endl;
                            buckets.push(new_bucket, v, [v](galois::InsertBag<std::uint64_t>& bag) { bag.push(v); });
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
                        //std::cout << "poo2\n";
                       // std::cout << "cur dist: " << cur_dist << "    new d " << new_dist << std::endl;

                        cur_dist = new_dist;
                        uint64_t new_bucket = cur_dist / DELTA;

                        //std::cout << "bucket: " << new_bucket << std::endl;
                        if (new_bucket != cur_bucket)
                            distributed_min.update(new_bucket);

                        if (new_bucket == cur_bucket && req_type == LIGHT) {
                            sbb.push(v);
                        } else {
                            buckets.push(new_bucket, v, [v](galois::InsertBag<std::uint64_t>& bag) { bag.push(v); });
                        }
                    }
                }
            }
    },
    galois::loopname("relax_edges"), galois::steal()
    );
    //std::cout << "made it out" << std::endl;
    return any_updated;
 }

 void delta_step(Graph& graph, std::uint64_t src, BucketContainer& buckets) {
    galois::GReduceMin<std::uint64_t> distributed_min; // starts with uint64_t max
    std::uint64_t min_bucket = 0; // find intl
    galois::InsertBag<std::uint64_t> dummy_bag;


   // std::cout << "START\n";

    relax_edges_for_vertex(graph,buckets,0, BOTH, distributed_min, min_bucket, dummy_bag);
  //  std::cout << "out\n";
    
    if (buckets[min_bucket].empty()) {
        min_bucket = distributed_min.reduce();
       // std::cout << "First min bucket " << min_bucket << std::endl;
    }

    distributed_min.reset();
    std::cout << "LESS GET STARTED\n";
    // only relax the minimum bucket
    while (true)  { // fix cond
        auto& bucket = buckets[min_bucket]; // retrieves min insert bag, ref
        std::unordered_set<std::uint64_t> r;
        
        bool any_updated_light;
        // std::cout << "min bucket: " << min_bucket << std::endl;
        
        while (!bucket.empty()) {
          //  std::cout << "NOT EMPTY FOR BUCKET " << min_bucket << std::endl;
            uint64_t cnt = 0;
            galois::InsertBag<std::uint64_t> cur_bucket_bag;
            for (std::uint64_t vertex : bucket) {
            //    std::cout << "VERTEX: " << vertex << std::endl;
                any_updated_light = relax_edges_for_vertex(graph,buckets, vertex,LIGHT, distributed_min, min_bucket, cur_bucket_bag);
                r.insert(vertex);
                cnt++;
            }
            bucket = std::move(cur_bucket_bag);
           // std::cout << "cnt " << cnt << std::endl;
        }

        // now the bucket is empty, we can free the assoc heap memory
        //std::cout << "erasing " << std::endl;
        buckets.erase(min_bucket);

       // std::cout << "ROUND 2\n";
        bool any_updated_heavy;
        for (std::uint64_t vertex : r) {
            any_updated_heavy = relax_edges_for_vertex(graph,buckets, vertex,HEAVY, distributed_min, min_bucket, dummy_bag); // will; now refill bucket 
        }
        if (!any_updated_heavy && !any_updated_light) {
            std::cout << "FINAL STEP\n";
            uint64_t test = ~0U;
            for (uint64_t i = 0; i < NUM_BUCKETS; i++) {
                if (!buckets.is_empty(i)) {
                    test = std::min(test, buckets.get_min(i));
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
