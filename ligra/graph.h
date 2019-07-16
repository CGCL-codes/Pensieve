#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "vertex.h"
// #include "compressedVertex.h"
#include "parallel.h"
#include "get_mem.h"

using namespace std;

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

// Class that handles implementation specific freeing of memory
// owned by the graph
struct Deletable {
public:
  virtual void del() = 0;
};

struct Empty_Mem : public Deletable {
  void del() {}
};

template <class vertex>
struct Uncompressed_Mem : public Deletable {
public:
  vector<vertex> V;
  long n;
  long m;
  void* allocatedInplace, * inEdges;

  Uncompressed_Mem(vector<vertex> VV, long nn, long mm, void* ai, void* _inEdges = NULL)
  : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges) { }
  Uncompressed_Mem(vertex* VV, long nn, long mm, void* ai, void* _inEdges = NULL)
  :  n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges) {
    for(int i=0; i<nn; i++) {
      V.push_back(*(VV+i));
    }
  }
  void del() {
    if (allocatedInplace == NULL)
      for (long i=0; i < n; i++) V[i].del();
    // else free(allocatedInplace);
    V.clear();
    if(inEdges != NULL) free(inEdges);
  }
};

template <class vertex>
struct Compressed_Mem : public Deletable {
public:
  vertex* V;
  char* s;

  Compressed_Mem(vertex* _V, char* _s) :
                 V(_V), s(_s) { }

  void del() {
    free(V);
    free(s);
  }
};

template <class vertex>
struct graph {
private:
  int max_version;
  vector<vertex> V;
  int version;//version id
public:
  long n;
  long m;
  bool transposed;
  bool _is_hybrid = false;
  uintE* flags;
  Deletable *D;

  graph(vertex* _V, long _n, long _m, Deletable* _D) : n(_n), m(_m),
    D(_D), flags(NULL), transposed(0), version(0) {
    for(int i=0; i<_n; i++) {
      V.push_back(*(_V+i));
    }
  }

  graph(vector<vertex>& _V, long _n, long _m, Deletable* _D) : V(_V), n(_n), m(_m),
  D(_D), flags(NULL), transposed(0), version(0) {}

  graph(vector<vertex>& _V, long _n, long _m, Deletable* _D, bool _hybrid) : V(_V), n(_n), m(_m),
  D(_D), flags(NULL), transposed(0), version(0), _is_hybrid(_hybrid) {}

  graph(vertex* _V, long _n, long _m, Deletable* _D, uintE* _flags, int _version) : V(_V),
  n(_n), m(_m), D(_D), flags(_flags), transposed(0), version(_version) {}

  ~graph() {}

  void del() {
    if (flags != NULL) free(flags);
    delete D;
  }

  void transpose() {
    // if ((typeid(vertex) == typeid(asymmetricVertex)) ||
    //     (typeid(vertex) == typeid(compressedAsymmetricVertex))) {
    //   parallel_for(long i=0;i<n;i++) {
    //     V[i].flipEdges();
    //   }
    //   transposed = !transposed;
    // }
  }
  vertex* getvertex() {
	  return V.data();
  }

  inline vertex * getvertex(uintT j) {
    if (j >= n) {
      resize(j+1);
    }
    return &V[j];
  }

  void resize(uintT new_size) {
    cout << "vertex expand from " << n << " to " << new_size << endl;
    vertex tmp;
    // use tmp to init other than possible random number.
    V.resize(new_size, tmp);
    n = V.size();
  }

  int get_edge_number() {
    return m;
  }

  int get_edge_capicity() {
    int count = 0;
    for (auto v : V) {
      count += v.outNeighbors.get_cap();
    }
    return count;
  }

  int get_version() {
    return version;
  }

  int get_max_version() {
    return max_version;
  }

  void update_m(void) {
    int ret = 0;
    for (auto v : V) {
      ret += v.getInDegree();
    }
    m = ret;
  }

  void add_m(int a) {
    if (m+a < 0) {
      update_m();
    }
    m += a;
  }

  void set_version(int vers) {
    version = vers;
    if (version > max_version) {
      max_version = version;
    }
  }

  uintE accessAllEdges() {
    uintE ret = 0;
    // int count = 0;
    // int length = 0;
    for (auto i = 0; i < n; i++) {
      ret ^= V[i].getDirtyData();
      // length = V[i].getInDegree();
      // vertex & tmp = V[i];
      // // // auto tmp = V[i].getInNeighbors();
      // for (auto j = 0; j < length; j++) {
      //   ret ^= tmp.getInNeighbor(j);
      //   // ret ^= *(tmp+j);
      // }
      // count += length;
    }
    // cout << "count " << n << " vertices and " << count << " edges " << endl;
    return ret;
  }

  uintE access_vertex(uintE vtx) {
    uintT ret = 0;
    for (auto i=0; i<V[vtx].getInDegree(); i++) {
      ret ^= V[vtx].getInNeighbor(i);
    }
    return ret;
  }

  bool is_hybrid() {
    return _is_hybrid;
  }

  bool is_high_degree_vertex(uintT index) {
    return V[index].is_high_degree_vertex();
  }

  void shrink() {
    for (auto i: V) {
      i.shrink();
    }
  }
};

#endif
