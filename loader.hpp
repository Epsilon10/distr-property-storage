#pragma once

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


#include "galois/Galois.h"
#include "edge_data.hpp"
#include "galois/graphs/Ben_CSR_64.h"

std::size_t constexpr NUM_NODES_LOG_2 = 20;
std::size_t constexpr NUM_NODES = 1 << NUM_NODES_LOG_2;
std::size_t constexpr EDGE_FACTOR = 16;
std::size_t constexpr NUM_EDGES = EDGE_FACTOR * NUM_NODES;
std::size_t constexpr NUM_WORKER_THREADS = 6; // TODO: increase
std::size_t num_merged = 0;
std::size_t constexpr DATA_SIZE_BYTES = 64;

using EdgeData_N = EdgeData_NBytes<DATA_SIZE_BYTES>;

struct NodeData {
    std::uint64_t dist;
    bool updated{false};
};
//using EdgeData_N = std::uint64_t;

// std::size_t constexpr NUM_EDGES_PER_THREAD = NUM_EDGES / NUM_WORKER_THREADS;

std::string const BASE_PATH = "/var/local/graphs/ben/tg_2/";
// std::string const BASE_PATH = "/var/local/graphs/ben/bg_2/";
//std::string const FILENAME_BASE = BASE_PATH + "graph_30_16";
std::string const FILENAME_BASE = BASE_PATH + "test";

using node_map_t = std::unordered_map<std::pair<std::uint64_t, std::uint64_t>, EdgeData_N, hash_pair>;

//node_map_t global_node_map;
std::unordered_map<std::pair<std::uint64_t, std::uint64_t>, std::size_t, hash_pair> global_node_map;
std::vector<EdgeData_N> edge_data_arr(NUM_EDGES);

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
    size_t w_no = 0;
    while (true) {
        node_map_t local_map;
        if (!wait_and_pop(local_map)) {
            // queue empty and terminate signal sent
            break;
        }
        
        for (auto [edge, weight] : local_map) {
            std::uint64_t src = edge.first;
            std::uint64_t dst = edge.second;
           //  std::cout << "WEIGHT THING " << weight << std::endl;
           
            //global_node_map[std::make_pair(src, dst)] = weight;
           global_node_map[std::make_pair(src, dst)] = w_no;
           edge_data_arr[w_no] = EdgeData_N();
            //std::cout << "WEIGHT NUMBER " << w_no << std::endl;


            edge_list[src].push_back(dst);
            w_no++;
        }
        
    }
}

void worker_thread_func(std::size_t thread_id) {
    node_map_t local_map;

    std::size_t const num_process_per_epoch = num_process_per_epoch_array[thread_id];


    std::string filename = FILENAME_BASE + "-" + std::to_string(thread_id) + ".el";
   // std::string filename = FILENAME_BASE + "-" + "0"+std::to_string(thread_id) + ".txt";
    std::cout << "FILENAME: " << filename << std::endl;
    
    std::string line;
    std::size_t e_no = 0;
    std::size_t e1, e2;
    char x;
    std::size_t e_no_this_epoch = 0;
    //std::size_t weight;

    std::ifstream infile(filename);

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        iss >> e1 >> x >> e2;
        

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

auto load_csr_64_graph() {    
    thread_vec.reserve(NUM_WORKER_THREADS);

    std::thread merger_thread(merger_thread_func);

    for (std::size_t i = 0; i < NUM_WORKER_THREADS; i++) {
        thread_vec.emplace_back(worker_thread_func, i);
    }
    std::cout << "STARTING COUNTDOWN\n";
    for (auto& thread : thread_vec) {
        thread.join();
        std::cout << "JOIN\n";
    }

    all_threads_done = true;

    std::cout << "COUNTDOWN OVER\n";
    work_queue_cv.notify_all();
    merger_thread.join();

    std::cout << "MAP SIZE: " << global_node_map.size() << std::endl;
 /*
    return galois::graphs::Ben_CSR_64<std::uint64_t, EdgeData_N, false>{
        NUM_NODES, NUM_EDGES, std::move(global_node_map), std::move(edge_list)
    };*/
    
    
    return galois::graphs::Ben_CSR_64<std::uint64_t, EdgeData_N, true>{
        NUM_NODES, NUM_EDGES, std::move(global_node_map),std::move(edge_data_arr), std::move(edge_list)
    };
}

std::uint64_t get_random_node() {
    int rand = std::rand();
    auto it = std::next(global_node_map.begin(), rand % global_node_map.size());
    std::uint64_t node = (it->first).first;
    return edge_list[node].size() >= 1 ? node : get_random_node();
}
