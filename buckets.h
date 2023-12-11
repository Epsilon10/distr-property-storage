#pragma once

#include <vector>
#include <mutex>
#include <cstdlib>
#include "queue.h"

struct BucketContainer {
    struct Vertex{
        uint64_t val;
        Vertex* next = nullptr;
        Vertex(uint64_t val) : val(val) { }
    };

    struct Bucket {
        uint64_t number{0};
        Queue<Vertex> vertices;
        Bucket* next = nullptr;
        Bucket(uint64_t number) : number(number) {}
    };

    uint64_t capacity;
    std::vector<Queue<Bucket>> table;
    std::vector<uint64_t> min_bucket_per_slot;
    std::vector<std::mutex> locks;

    BucketContainer(uint64_t capacity) : 
        capacity(capacity), table(capacity), min_bucket_per_slot(capacity, ~0U), locks(capacity) { }

    // never called inside parallel section
    Bucket* get(uint64_t bucket_no) {
        std::uint64_t idx = bucket_no & (capacity - 1);
        auto& slot_q = table[idx];
        
        Bucket* bucket = slot_q.get([bucket_no](Bucket* bucket) { return bucket->number == bucket_no; });

        assert (bucket != nullptr);
        
        return bucket;
    }

    // can be called in parallel section
    void push(uint64_t bucket_no, uint64_t vertex) {
        std::uint64_t idx = bucket_no & (capacity - 1);
     //   std::cout << "Acquring lock " << idx << " for indice " << i << std::endl;
        std::lock_guard g{locks[idx]};

        auto& slot_q = table[idx];
        Bucket* bucket = slot_q.get([bucket_no](Bucket* bucket) { return bucket->number == bucket_no; });

        if (bucket != nullptr) {
            bucket->vertices.add(new Vertex(vertex));
        } else {
            bucket = new Bucket(bucket_no);
            bucket->vertices.add(new Vertex(vertex));
            slot_q.add(bucket);
            min_bucket_per_slot[idx] = std::min(min_bucket_per_slot[idx], bucket_no);
        }

    }

    // called in sync substrate which is per    
    void remove_vertex_from_bucket(uint64_t bucket_no, uint64_t vertex) {
        std::uint64_t idx = bucket_no & (capacity - 1);
     //   std::cout << "Acquring lock " << idx << " for indice " << i << std::endl;
        std::lock_guard g{locks[idx]};

        auto& slot_q = table[idx];
        Bucket* bucket = slot_q.get([bucket_no](Bucket* bucket) { return bucket->number == bucket_no; });
        assert(bucket != nullptr);

        bucket->vertices.remove_first([vertex](Vertex* vtx) { return vtx->val == vertex; });
    }

    // never called in parallel section
    void remove_bucket(uint64_t bucket_no) {
        std::uint64_t idx = bucket_no & (capacity - 1);
        auto& slot_q = table[idx];

        slot_q.remove_first([bucket_no](Bucket* bucket) { return bucket->number == bucket_no; });

        update_slot_min_bucket(bucket_no);
    }

    // never called inside parallel section
    bool is_empty(uint64_t bucket_no) {
        std::uint64_t idx = bucket_no & (capacity - 1);
        auto& slot_q = table[idx];

        Bucket* bucket = slot_q.get([bucket_no](Bucket* bucket) { return bucket->number == bucket_no; });
        return bucket == nullptr || bucket->vertices.is_empty();
    }

    void update_slot_min_bucket(uint64_t idx) {
        auto& slot_q = table[idx];

        uint64_t min_buck = ~0U;
        slot_q.all([&min_buck](Bucket* bucket) { min_buck = std::min(min_buck, bucket->number); });
        min_bucket_per_slot[idx] = min_buck;
    }

    uint64_t get_min(uint64_t idx) {
        return min_bucket_per_slot[idx];
    }
};
