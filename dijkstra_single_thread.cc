#include <iostream>
#include <vector>
#include <thread>
#include <cstddef>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <queue>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <chrono>
#include <iterator>
#include <cstdlib>
#include <libmemcached/memcached.h>


#include "galois/Galois.h"
#include "edge_data.hpp"
#include "galois/graphs/Ben_CSR_64.h"

std::size_t constexpr NUM_NODES_LOG_2 = 16;
std::size_t constexpr NUM_NODES = 1 << NUM_NODES_LOG_2;
std::size_t constexpr EDGE_FACTOR = 16;
std::size_t constexpr NUM_EDGES = EDGE_FACTOR * NUM_NODES;
std::size_t constexpr NUM_WORKER_THREADS = 4; // TODO: increase
std::size_t num_merged = 0;
std::size_t constexpr DATA_SIZE_BYTES = 8;

using EdgeData_N = EdgeData_NBytes<DATA_SIZE_BYTES>;


// std::size_t constexpr NUM_EDGES_PER_THREAD = NUM_EDGES / NUM_WORKER_THREADS;

std::string const BASE_PATH = "/var/local/graphs/ben/";
//std::string const FILENAME_BASE = BASE_PATH + "graph_30_16";
std::string const FILENAME_BASE = BASE_PATH + "output16";

using node_map_t = std::unordered_map<std::pair<std::uint64_t, std::uint64_t>, EdgeData_N, hash_pair>;

node_map_t global_node_map;

std::array<uint64_t, NUM_WORKER_THREADS> num_process_per_epoch_array = {{10000, 20000, 30000, 40000}};

std::vector<std::vector<std::uint64_t>> edge_list(NUM_NODES);

std::vector<std::thread> thread_vec;

// thread safe work queue variables
std::queue<node_map_t> work_queue;
std::mutex work_queue_mutex;
std::condition_variable work_queue_cv;

bool all_threads_done = false;

std::uint64_t get_random_node();
void my2_memcached_test();

void work_queue_push(node_map_t&& local_map) {
    std::lock_guard<std::mutex> lock(work_queue_mutex);
    work_queue.push(std::move(local_map));
    work_queue_cv.notify_one();
}

bool work_queue_is_empty() {
    std::lock_guard<std::mutex> lock(work_queue_mutex);
    return work_queue.empty();
}

bool wait_and_pop(node_map_t& local_map) {
    std::unique_lock<std::mutex> lock(work_queue_mutex);
    work_queue_cv.wait(lock, []{ return !work_queue.empty() || all_threads_done; });
    if (work_queue.empty()) { return false; }

    local_map = std::move(work_queue.front());
    work_queue.pop();
    lock.unlock();    
    return true;
}

void merger_thread_func() {
    while (true) {
        node_map_t local_map;
        if (!wait_and_pop(local_map)) {
            // queue empty and terminate signal sent
            break;
        }
        
        for (auto [edge, weight] : local_map) {
            std::uint64_t src = edge.first;
            std::uint64_t dst = edge.second;
            
            global_node_map[std::make_pair(src, dst)] = weight;
            edge_list[src].push_back(dst);
        }
        
    }
}

void worker_thread_func(std::size_t thread_id) {
    node_map_t local_map;

    std::size_t const num_process_per_epoch = num_process_per_epoch_array[thread_id];


    //std::string filename = FILENAME_BASE + "-" + std::to_string(thread_id) + ".el";
    std::string filename = FILENAME_BASE + "-" + "0"+std::to_string(thread_id) + ".txt";
    std::cout << "FILENAME: " << filename << std::endl;
    
    std::string line;
    std::size_t e_no = 0;
    std::size_t e1, e2;
    std::size_t e_no_this_epoch = 0;
    //std::size_t weight;

    std::ifstream infile(filename);

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        iss >> e1 >> e2;

        // e2 as dst edge, e_no is edge weight
        local_map[std::make_pair(e1, e2)] = EdgeData_N();

        e_no++;
        e_no_this_epoch++;

        if (e_no_this_epoch == num_process_per_epoch) {
            work_queue_push(std::move(local_map));
            local_map.clear();
           //std::cout << "LMAP SIZE: " << local_map.size() << std::endl;
            e_no_this_epoch = 0;
        }        
    }

    work_queue_push(std::move(local_map));
    local_map.clear();

    std::cout << "E_NO: " << e_no << std::endl; 
}

int main() {
    galois::SharedMemSys G;
    my2_memcached_test();
    std::cout << "NUM EDGES: " << NUM_EDGES << std::endl; 
    memcached_st* m = memcached_create(NULL);
    std::cout << "NEW\n";
    
    thread_vec.reserve(NUM_WORKER_THREADS);

    std::thread merger_thread(merger_thread_func);

    for (std::size_t i = 0; i < NUM_WORKER_THREADS; i++) {
        thread_vec.emplace_back(worker_thread_func, i);
    }
    std::cout << "STARTING COUNTDOWN\n";
    for (auto& thread : thread_vec) {
        thread.join();
    }

    all_threads_done = true;

    std::cout << "COUNTDOWN OVER\n";
    work_queue_cv.notify_all();
    merger_thread.join();

    std::cout << "MAP SIZE: " << global_node_map.size() << std::endl;
    
    /*
    for (auto [key, val] : global_node_map) {
        std::cout << key.first << " " << key.second << " " << val.weight << std::endl;
    }   */

    std::array<std::uint64_t, 64> random_search_keys;

    for (std::size_t i = 0; i < 64; i++) {
        random_search_keys[i] = get_random_node();
    }
    
    galois::graphs::Ben_CSR_64<std::uint64_t, EdgeData_N, false> graph{NUM_NODES, NUM_EDGES, std::move(global_node_map)
        ,std::move(edge_list)};
    
    global_node_map.clear();

    std::cout << "BEGINNING KERNEL 2: " << std::endl;
    float avg_time = 0.0;
    for (std::size_t i = 0; i < 64; i++) {
        std::vector<std::size_t> distance(NUM_NODES);
        for (std::size_t i = 0; i < NUM_NODES; i++) {
            distance[i] = std::numeric_limits<std::size_t>::max();
        }
        std::uint64_t start_node = random_search_keys[i];
        std::cout << "RANDOM NODE SELECTED: " << start_node << std::endl;
        distance[start_node] = 0;

        auto start = std::chrono::high_resolution_clock::now();

        std::priority_queue<std::pair<std::size_t, std::size_t>> q;

        q.push(std::make_pair(0, start_node));
        std::vector<bool> processed(NUM_NODES, false);

        while (!q.empty()) {
            std::uint64_t a = q.top().second;
            q.pop();
            if (processed[a]) continue;

            processed[a] = true;

            auto degree = graph.getDegree(a);
            auto edge_itr_start = graph.edge_begin(a);
            auto edge_itr_end = graph.edge_end(a);

            for (std::size_t i = 0; i < degree; i++) {
                auto b = graph.getEdgeDst(edge_itr_start + i);
                auto& w = graph.getEdgeData(a, b);

                if (w +distance[a] < distance[b]) {
                    distance[b] = w+distance[a];
                    q.push(std::make_pair(-distance[b],b));
                }
            }
            
        }

        auto stop = std::chrono::high_resolution_clock::now();
        auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    
    /*
        for (std::size_t i = 1; i < NUM_NODES; i++) {
            std::cout << "DISTANCE FROM SRC: " << distance[i] << std::endl;
        }*/
        
        std::cout << "TRIAL: " << i << " exec time: " << time_ms.count() << std::endl;
        avg_time += time_ms.count();
    }
    std::cout << "MEAN TIME: " << avg_time / 64 << std::endl;
}

std::uint64_t get_random_node() {
    int rand = std::rand();
    auto it = std::next(global_node_map.begin(), rand % global_node_map.size());
    std::uint64_t node = (it->first).first;
    return edge_list[node].size() >= 1 ? node : get_random_node();
}

void my2_memcached_test() {
  memcached_server_st *servers = NULL;
  memcached_st *memc;
  memcached_return rc;
  char *key = "keystring";
  char *value = "keyvalue";

  char *retrieved_value;
  size_t value_length;
  uint32_t flags;

  memc = memcached_create(NULL);
  servers = memcached_server_list_append(servers, "localhost", 11211, &rc);
  rc = memcached_server_push(memc, servers);

  if (rc == MEMCACHED_SUCCESS)
    fprintf(stderr, "Added server successfully\n");
  else
    fprintf(stderr, "Couldn't add server: %s\n", memcached_strerror(memc, rc));

  rc = memcached_set(memc, key, strlen(key), value, strlen(value), (time_t)0, (uint32_t)0);

  if (rc == MEMCACHED_SUCCESS)
    fprintf(stderr, "Key stored successfully\n");
  else
    fprintf(stderr, "Couldn't store key: %s\n", memcached_strerror(memc, rc));

  retrieved_value = memcached_get(memc, key, strlen(key), &value_length, &flags, &rc);
  printf("Yay!\n");

  if (rc == MEMCACHED_SUCCESS) {
    fprintf(stderr, "Key retrieved successfully\n");
    printf("The key '%s' returned value '%s'.\n", key, retrieved_value);
    free(retrieved_value);
  }
  else
    fprintf(stderr, "Couldn't retrieve key: %s\n", memcached_strerror(memc, rc));

}