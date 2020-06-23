#ifndef VERTEX_H
#define VERTEX_H

#include <chrono>
#include <algorithm>
#include "vertexSubset.h"
#include <vector>
#include "myutil.h"
#include "copyed_array.h"
using namespace std;

namespace decode_uncompressed {

  // Used by edgeMapDense. Callers ensure cond(v_id). For each vertex, decode
  // its in-edges, and check to see whether this neighbor is in the current
  // frontier, calling update if it is. If processing the edges sequentially,
  // break once !cond(v_id).
  template <class vertex, class F, class G, class VS>
  inline void decodeInNghBreakEarly(vertex* v, long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    uintE d = v->getInDegree();
    if (!parallel || d < 1000) {
      for (size_t j=0; j<d; j++) {
        uintE ngh = v->getInNeighbor(j);
        if (vertexSubset.isIn(ngh)) {
#ifndef WEIGHTED
          auto m = f.update(ngh, v_id);
#else
          auto m = f.update(ngh, v_id, v->getInWeight(j));
#endif
          g(v_id, m);
        }
        if(!f.cond(v_id)) break;
      }
    } else {
      parallel_for(size_t j=0; j<d; j++) {
        uintE ngh = v->getInNeighbor(j);
        if (vertexSubset.isIn(ngh)) {
#ifndef WEIGHTED
          auto m = f.updateAtomic(ngh, v_id);
#else
          auto m = f.updateAtomic(ngh, v_id, v->getInWeight(j));
#endif
          g(v_id, m);
        }
      }
    }
  }

  // Used by edgeMapDenseForward. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <class V, class F, class G>
  inline void decodeOutNgh(V* v, long i, F &f, G &g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
#ifndef WEIGHTED
      auto m = f.updateAtomic(i,ngh);
#else
      auto m = f.updateAtomic(i,ngh,v->getOutWeight(j));
#endif
        g(ngh, m);
      }
    });
  }

  // Used by edgeMapSparse. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <class V, class F, class G>
  inline void decodeOutNghSparse(V* v, long i, uintT o, F &f, G &g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
#ifndef WEIGHTED
        auto m = f.updateAtomic(i, ngh);
#else
        auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
        g(ngh, o+j, m);
      } else {
        g(ngh, o+j);
      }
    });
  }

  // Used by edgeMapSparse_no_filter. Sequentially decode the out-neighbors,
  // and compactly write all neighbors satisfying g().
  template <class V, class F, class G>
  inline size_t decodeOutNghSparseSeq(V* v, long i, uintT o, F &f, G &g) {
    uintE d = v->getOutDegree();
    size_t k = 0;
    for (size_t j=0; j<d; j++) {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
#ifndef WEIGHTED
        auto m = f.updateAtomic(i, ngh);
#else
        auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
        bool wrote = g(ngh, o+k, m);
        if (wrote) { k++; }
      }
    }
    return k;
  }

  // Decode the out-neighbors of v, and return the number of neighbors
  // that satisfy f.
  template <class V, class F>
  inline size_t countOutNgh(V* v, long vtx_id, F& f) {
    uintE d = v->getOutDegree();
    if (d < 2000) {
      size_t ct = 0;
      for (size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
#ifndef WEIGHTED
        if (f(vtx_id, ngh))
#else
        if (f(vtx_id, ngh, v->getOutWeight(i)))
#endif
          ct++;
      }
      return ct;
    } else {
      size_t b_size = 2000;
      size_t blocks = 1 + ((d-1)/b_size);
      auto cts = array_imap<uintE>(blocks, [&] (size_t i) { return 0; });
      parallel_for_1(size_t i=0; i<blocks; i++) {
        size_t s = b_size*i;
        size_t e = std::min(s + b_size, (size_t)d);
        uintE ct = 0;
        for (size_t j = s; j < e; j++) {
          uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
          if (f(vtx_id, ngh))
#else
          if (f(vtx_id, ngh, v->getOutNeighbor(j)))
#endif
            ct++;
        }
        cts[i] = ct;
      }
      size_t count = 0;
      return pbbs::reduce_add(cts);
    }
  }

  // Decode the out-neighbors of v. Apply f(src, ngh) and store the result
  // using g.
  template <class V, class E, class F, class G>
  inline void copyOutNgh(V* v, long src, uintT o, F& f, G& g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
#ifdef WEIGHTED
      E val = f(src, ngh, v->getOutWeight(j));
#else
      E val = f(src, ngh);
#endif
      g(ngh, o+j, val);
    });
  }

  // TODO(laxmand): Add support for weighted graphs.
  template <class V, class Pred>
  inline size_t packOutNgh(V* v, long vtx_id, Pred& p, bool* bits, uintE* tmp) {
    uintE d = v->getOutDegree();
    if (d < 5000) {
      size_t k = 0;
      for (size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
        if (p(vtx_id, ngh)) {
          v->setOutNeighbor(k, ngh);
          k++;
        }
      }
      v->setOutDegree(k);
      return k;
    } else {
      parallel_for(size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
        tmp[i] = ngh;
        bits[i] = p(vtx_id, ngh);
      }
      size_t k = sequence::pack(tmp, v->getOutNeighbors(), bits, d);
      v->setOutDegree(k);
      return k;
    }
  }

}

struct HVertex
{
  // when store edges of the vertex
  vector<uintE> outNeighbors;

  // pos points to offset in carray; if pos == -1, it points to outNeighbors.
  int pos = -1;
  bool _is_high_degree_vertex = true;

  static copyed_array carray;

  // init func when a graph append more vertex during a few version later.
  // give these vertex a empty version 0 to make it valid
  HVertex() {
    pos = -1;
    _is_high_degree_vertex = true;
  }

  // if give 
  HVertex(vector<uintE>& n, bool _high_degree_vertex = true) {
    outNeighbors = n;
    _is_high_degree_vertex = _high_degree_vertex;
  }

  HVertex(uintE *data, uintT length, bool _high_degree_vertex=true) {
    outNeighbors.insert(outNeighbors.end(), data, data + length);
    _is_high_degree_vertex = _high_degree_vertex;
  }

  void setInNeighbors(uintE* _i) { return ; }
  void setOutNeighbors(uintE* _i) { return ; }
  void setInDegree(uintT _d) {  }
  void setOutDegree(uintT _d) { setInDegree(_d); }
  void flipEdges() {}

  uintE* getInNeighbors() {
    if (_is_high_degree_vertex || pos < 0) {
      // return outNeighbors.data() + pointers[pos];
      return outNeighbors.data();
    } else {
      return carray.edges.data() + carray.offsets[pos];
    }
  }
  uintE * getOutNeighbors() { return getInNeighbors(); }

  vector<uintE>  getLastVersion() {
    if (pos == -1) return outNeighbors;
    return vector<uintE> (carray.edges.begin() + carray.offsets[pos], carray.edges.begin() + carray.offsets[pos+1]);
  }

  void index_delete(uintT index) {
    if (_is_high_degree_vertex) {
      vector_index_delete(outNeighbors, index);
    }
  }

  void index_addtion(uintE data, uintT index) {
    if (_is_high_degree_vertex) {
      vector_index_addition(outNeighbors, index, data);
    }
  }

  void push_back(uintE data) {
    if (_is_high_degree_vertex) {
      outNeighbors.push_back(data);
    }
  }

  void push_back(uintE * data, uintT length) {
    if (_is_high_degree_vertex) {
      outNeighbors.insert(outNeighbors.end(), data, data+length);
    }
  }

  void push_back(vector<uintE>::iterator start, vector<uintE>::iterator end) {
    if (_is_high_degree_vertex) {
      outNeighbors.insert(outNeighbors.end(), start, end);
    }
  }

  void pop_back() {
    if (_is_high_degree_vertex) {
      outNeighbors.pop_back();
    }
  }

  bool is_high_degree_vertex() {
    return _is_high_degree_vertex;
  }

  void shrink() {
    if (is_high_degree_vertex()) {
      outNeighbors.shrink_to_fit();
    }
  }

  void sort() {
    if (_is_high_degree_vertex) {
      std::sort(outNeighbors.begin(), outNeighbors.end());
    }
  }

  // for low degree vertex to append a new version.
  void append(vector<uintE> &data, uintT vtx, uintT version) {
    if (_is_high_degree_vertex) {
      cout << "append to high vertex " << endl;
      return ;
    }
    int next_pos = carray.append(data, vtx, version, pos);
    
    if (pos < 0) {
      pos = -1*next_pos-1;
    }
  }

  void empty_step(uintT arg) {}

  void empty_step() {}

  void forward() {
    if (!_is_high_degree_vertex) {
      if (pos < 0) {
        pos = -1*(pos+1);
      } else {
        pos = carray.nexts[pos];
      }
    }
  }

  // void forward(uintT ver) {
  //   if (carray.versions[])
  // }

  void backward() {
    if (!_is_high_degree_vertex) {
      pos = carray.froms[pos];
    }
  }

  void reserve(uintT s) {
    if (_is_high_degree_vertex) {
      outNeighbors.reserve(s);
    }
  }

  void setInNeighbor(uintT j, uintE ngh) {}
  void setOutNeighbor(uintT j, uintE ngh) {}
  uintE getInNeighbor(uintT j) {
    if (_is_high_degree_vertex || pos < 0) {
      return outNeighbors[j];
    } else {
      return carray.edges[carray.offsets[pos]+j];
    }
  }
  uintE getOutNeighbor(uintT j) {
    return getInNeighbor(j);
  }

  uintT getInDegree() {
    if (_is_high_degree_vertex || pos < 0) {
      return outNeighbors.size();
    } else {
      return carray.offsets[pos+1] - carray.offsets[pos];
    }
  }
  uintT getOutDegree() {
    return getInDegree();
  }

  uintT getDirtyData() {
    uintT ret = 0;
    auto start = getInNeighbors();
    for (auto i=0; i<getInDegree(); i++) {
      ret ^= *(start+i);
    }
    return ret;
  }

  void print() {
    cout << "print begin " << endl;
    for (auto i=0; i<getInDegree(); i++) {
      cout << *(getInNeighbors() + i) << " ";
    }
    cout << endl;
    cout << "print end " << endl;
  }

  int	find_val(const uintE &val) {
    auto start = getInNeighbors();
    auto end = getInNeighbors() + getInDegree();
    auto it = find(start, end, val);
    if (it == end) {
      return -1;
    }
    return it - start;
  }

  HVertex operator=(const HVertex & other) {
    if (&other != this) {
      outNeighbors = other.outNeighbors;
      _is_high_degree_vertex = other._is_high_degree_vertex;
      pos = other.pos;
    }
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<HVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G& g) {
     decode_uncompressed::decodeOutNgh<HVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_uncompressed::decodeOutNghSparse<HVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_uncompressed::decodeOutNghSparseSeq<HVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_uncompressed::copyOutNgh<HVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_uncompressed::countOutNgh<HVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_uncompressed::packOutNgh<HVertex, F>(this, i, f, bits, tmp1);
  }
};

copyed_array HVertex::carray = copyed_array();

struct symmetricVertex {
#ifndef WEIGHTED
  vector<uintE> outNeighbors;
#else
  vector<intE> outNeighbors;
#endif
  void del() {}
#ifndef WEIGHTED
symmetricVertex(vector<uintE> n)
#else
symmetricVertex(vector<intE> n)
#endif
: outNeighbors(n) {}

symmetricVertex()
#ifndef WEIGHTED
{ outNeighbors = vector<uintE>(); }
#else
{ outNeighbors = vector<intE>(); }
#endif

#ifndef WEIGHTED
  uintE* getInNeighbors () { return outNeighbors.data(); }
  const uintE* getInNeighbors () const { return outNeighbors.data(); }
  uintE* getOutNeighbors () { return outNeighbors.data(); }
  const uintE* getOutNeighbors () const { return outNeighbors.data(); }
  uintE getInNeighbor(uintT j) const { return outNeighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }

  void setInNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  void setInNeighbors(uintE* _i) { return ; }
  void setOutNeighbors(uintE* _i) { return ; }
#else
  //weights are stored in the entry after the neighbor ID
  //so size of neighbor list is twice the degree
  intE* getInNeighbors () { return outNeighbors.data(); }
  const intE* getInNeighbors () const { return outNeighbors.data(); }
  intE* getOutNeighbors () { return outNeighbors.data(); }
  const intE* getOutNeighbors () const { return outNeighbors.data(); }
  intE getInNeighbor(intT j) const { return outNeighbors[2*j]; }
  intE getOutNeighbor(intT j) const { return outNeighbors[2*j]; }
  intE getInWeight(intT j) const { return outNeighbors[2*j+1]; }
  intE getOutWeight(intT j) const { return outNeighbors[2*j+1]; }
  void setInNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setInWeight(uintT j, intE wgh) { outNeighbors[2*j+1] = wgh; }
  void setOutWeight(uintT j, intE wgh) { outNeighbors[2*j+1] = wgh; }
  void setInNeighbors(intE* _i) {  }
  void setOutNeighbors(intE* _i) {  }
#endif

  void index_delete(uintT pos) { vector_index_delete(outNeighbors, pos); }
  void index_addtion(uintE data, uintT pos) { vector_index_addition(outNeighbors, pos, data); }
  void push_back(uintE data) { outNeighbors.push_back(data); }
  void push_back(vector<uintE>::iterator start, vector<uintE>::iterator end) { outNeighbors.insert(outNeighbors.end(), start, end); }
  void pop_back() { outNeighbors.pop_back(); }
  bool is_high_degree_vertex() { return true; }
  vector<uintE> getLastVersion() { return outNeighbors; }
  void forward(uintT ver) {}
  void forward() {}
  void backward(uintT ver) {}
  void backward() {}
  void reserve(uintT s) { outNeighbors.reserve(s); }
  void append(vector<uintE> data, uintT vtx, uintT version) {}
  void print() {
    cout << "print begin" << endl;
    for (auto i=0; i<getInDegree(); i++) {
      cout << getInNeighbor(i) << endl;
    }
    cout << "print end" << endl;
  }

  uintE getDirtyData() {
    uintE ret;
    for (auto i:outNeighbors) {
      ret ^= i;
    }
    return ret;
  }

  uintT getInDegree() const { return outNeighbors.size(); }
  uintT getOutDegree() const { return outNeighbors.size(); }
  void setInDegree(uintT _d) { return;  }
  void setOutDegree(uintT _d) { return;  }
  int	find_val(const uintE &val) { 
    auto it = find(outNeighbors.begin(), outNeighbors.end(), val);
    if (it == outNeighbors.end()) {
      return -1;
    } else {
      return it - outNeighbors.begin();
    }
  }
  void flipEdges() {}

  symmetricVertex& operator=(const symmetricVertex & other) {
    
    if (&other != this) {
      outNeighbors = other.outNeighbors;
    }
    return *this;
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<symmetricVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G& g) {
     decode_uncompressed::decodeOutNgh<symmetricVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_uncompressed::decodeOutNghSparse<symmetricVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_uncompressed::decodeOutNghSparseSeq<symmetricVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_uncompressed::copyOutNgh<symmetricVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_uncompressed::countOutNgh<symmetricVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_uncompressed::packOutNgh<symmetricVertex, F>(this, i, f, bits, tmp1);
  }
};

struct asymmetricVertex {
#ifndef WEIGHTED
  vector<uintE> inNeighbors, outNeighbors;
#else
  vector<intE> inNeighbors, outNeighbors;
#endif
  // uintT outDegree;
  // uintT inDegree;
  void del() {inNeighbors.clear(); outNeighbors.clear();}

#ifndef WEIGHTED
asymmetricVertex(vector<uintE> iN, vector<uintE> oN)
#else
asymmetricVertex(vector<intE> iN, vector<intE> oN)
#endif
: inNeighbors(iN), outNeighbors(oN){}

asymmetricVertex()
#ifndef WEIGHTED
{
	inNeighbors = vector<uintE>();
	outNeighbors = vector<uintE>();
}
#else
{	
	inNeighbors = vector<intE>();
	outNeighbors = vector<intE>();
}
#endif

#ifndef WEIGHTED
  uintE* getInNeighbors () { return inNeighbors.data(); }
  const uintE* getInNeighbors () const { return inNeighbors.data(); }
  uintE* getOutNeighbors () { return outNeighbors.data(); }
  const uintE* getOutNeighbors () const { return outNeighbors.data(); }
  uintE getInNeighbor(uintT j) const { return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  // TODO:: 
  void setInNeighbors(uintE* _i) { return;  }
  void setOutNeighbors(uintE* _i) { return;  }
#else
  intE* getInNeighbors () { return inNeighbors.data(); }
  const intE* getInNeighbors () const { return inNeighbors.data(); }
  intE* getOutNeighbors () { return outNeighbors.data(); }
  const intE* getOutNeighbors () const { return outNeighbors.data(); }
  intE getInNeighbor(uintT j) const { return inNeighbors[2*j]; }
  intE getOutNeighbor(uintT j) const { return outNeighbors[2*j]; }
  intE getInWeight(uintT j) const { return inNeighbors[2*j+1]; }
  intE getOutWeight(uintT j) const { return outNeighbors[2*j+1]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[2*j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setInWeight(uintT j, uintE wgh) { inNeighbors[2*j+1] = wgh; }
  void setOutWeight(uintT j, uintE wgh) { outNeighbors[2*j+1] = wgh; }
  void setInNeighbors(intE* _i) {  }
  void setOutNeighbors(intE* _i) {  }
#endif

  uintT getInDegree() const { return inNeighbors.size(); }
  uintT getOutDegree() const { return outNeighbors.size(); }
  uintT	find_val(const uintE &val) { return find(outNeighbors.begin(), outNeighbors.end(), val) - outNeighbors.begin(); }
  // TODO::::
  void setInDegree(uintT _d) { return;  }
  void setOutDegree(uintT _d) { return;  }
  void flipEdges() { inNeighbors.swap(outNeighbors);  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<asymmetricVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G &g) {
    decode_uncompressed::decodeOutNgh<asymmetricVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_uncompressed::decodeOutNghSparse<asymmetricVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_uncompressed::decodeOutNghSparseSeq<asymmetricVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_uncompressed::copyOutNgh<asymmetricVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_uncompressed::countOutNgh<asymmetricVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_uncompressed::packOutNgh<asymmetricVertex, F>(this, i, f, bits, tmp1);
  }

};

#endif
