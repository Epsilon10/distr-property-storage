#pragma once

#include <mutex>

template <typename T>
class Queue {
   public:
    T* volatile first = nullptr;
    T* volatile last = nullptr;
    
    Queue() : first(nullptr), last(nullptr) {}
    Queue(const Queue&) = delete;
    Queue& operator=(Queue&) = delete;
    Queue& operator=(Queue&& other) {
        this->first = other.first;
        this->last = other.last;

        other.first = nullptr;
        other.last = nullptr;
    }

    void add(T* t) {
        // g{lock};
        t->next = nullptr;
        if (first == nullptr) {
            first = t;
        } else {
            last->next = t;
        }
        last = t;
    }

    void add_stack(T* t) {
        // g{lock};
        t->next = first;
        first = t;
        if (last == nullptr) {
            last = t;
        }
    }

    void remove_all() {
        // g{lock};
         while (first != nullptr) {
            T* old = first;
            first = first->next;

            delete old;
        }

        first = nullptr;
        last = nullptr;
    }
    

    bool is_empty() {
        // g{lock};
        return first == nullptr;
    }

    template<typename Compare>
    void remove_first(Compare cmp) {
        // cmp = node->tag == i
        // g{lock};
        assert(first != nullptr);

        if (first->next == nullptr) {
            assert(cmp(first));
            delete first;
            first = nullptr;            
            return;
        }

        T* last = nullptr;
        T* node = first;
        while (node != nullptr) {
            if (cmp(node)) {
                last->next = node->next;
                delete node;
                break;
            }
            last = node;
            node = node->next;
        }
    }

    template<typename Compare>
    T* get(Compare cmp) {
        // g{lock};

        T* cur = first;
        if (cur == nullptr) return nullptr;

        while (cur != nullptr) {
            if (cmp(cur)) {
                return cur;
            }

            cur = cur->next;
        }

        return nullptr;
    }

    template<typename Work> 
    void all(Work w) {
        T* cur = first;
        
        while (cur != nullptr) {
            work(cur);
            cur = cur->next;
        }
    }
};

