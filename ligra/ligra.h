// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "vertex.h"
// #include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"
#include "index_map.h"
#include "edgeMap_utils.h"
#include "delta.h"
#include "exp_distribution.h"
#include "get_mem.h"
#include "big_delta.h"
#include "version_graph.h"
#include <chrono>
#include <time.h>
#include <vector>
using namespace std;

//*****START FRAMEWORK*****

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_no_filter = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
inline bool should_output(const flags& fl) { return !(fl & no_output); }

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDense(graph<vertex> * GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA->n;
  ChronoTimer cter;
  // cout << cter.elapsed_milli() << endl;
  // cout << should_output(fl) << endl;
  // vertex *G = GA->getvertex();
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_gen<data>(next);
    parallel_for (long v=0; v<n; v++) {
      std::get<0>(next[v]) = 0;
      if (f.cond(v)) {
        GA->getvertex(v)->decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
    // cout << cter.elapsed_milli() << endl;
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_nooutput_gen<data>();
    parallel_for (long v=0; v<n; v++) {
      if (f.cond(v)) {
        GA->getvertex(v)->decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDenseForward(graph<vertex> * GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA->n;
  // vertex *G = GA->getvertex();
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_forward_gen<data>(next);
    parallel_for(long i=0;i<n;i++) { std::get<0>(next[i]) = 0; }
    parallel_for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        GA->getvertex(i)->decodeOutNgh(i, f, g);
      }
    }
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_forward_nooutput_gen<data>();
    parallel_for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        GA->getvertex(i)->decodeOutNgh(i, f, g);
      }
    }
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse(graph<vertex>& GA, vertex* frontierVertices, VS& indices,
        uintT* degrees, uintT m, F &f, const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  S* outEdges;
  long outEdgeCount = 0;

  if (should_output(fl)) {
    uintT* offsets = degrees;
    outEdgeCount = sequence::plusScan(offsets, offsets, m);
    outEdges = newA(S, outEdgeCount);
    auto g = get_emsparse_gen<data>(outEdges);
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i), o = offsets[i];
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, o, f, g);
    }
  } else {
    auto g = get_emsparse_nooutput_gen<data>();
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i);
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, 0, f, g);
    }
  }

  if (should_output(fl)) {
    S* nextIndices = newA(S, outEdgeCount);
    if (fl & remove_duplicates) {
      if (GA.flags == NULL) {
        GA.flags = newA(uintE, n);
        parallel_for(long i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
      }
      auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(outEdges[i]); };
      remDuplicates(get_key, GA.flags, outEdgeCount, n);
    }
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
    free(outEdges);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  } else {
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse(graph<vertex>& GA, vertex** frontierVertices, VS& indices,
        uintT* degrees, uintT m, F &f, const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  S* outEdges;
  long outEdgeCount = 0;

  if (should_output(fl)) {
    uintT* offsets = degrees;
    outEdgeCount = sequence::plusScan(offsets, offsets, m);
    outEdges = newA(S, outEdgeCount);
    auto g = get_emsparse_gen<data>(outEdges);
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i), o = offsets[i];
      vertex* vert = frontierVertices[i];
      vert->decodeOutNghSparse(v, o, f, g);
    }
  } else {
    auto g = get_emsparse_nooutput_gen<data>();
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i);
      vertex* vert = frontierVertices[i];
      vert->decodeOutNghSparse(v, 0, f, g);
    }
  }

  if (should_output(fl)) {
    S* nextIndices = newA(S, outEdgeCount);
    if (fl & remove_duplicates) {
      if (GA.flags == NULL) {
        GA.flags = newA(uintE, n);
        parallel_for(long i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
      }
      auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(outEdges[i]); };
      remDuplicates(get_key, GA.flags, outEdgeCount, n);
    }
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
    free(outEdges);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  } else {
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter(graph<vertex>& GA,
    vertex* frontierVertices, VS& indices, uintT* offsets, uintT m, F& f,
    const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  long outEdgeCount = sequence::plusScan(offsets, offsets, m);
  S* outEdges = newA(S, outEdgeCount);

  auto g = get_emsparse_no_filter_gen<data>(outEdges);

  // binary-search into scan to map workers->chunks
  size_t b_size = 10000;
  size_t n_blocks = nblocks(outEdgeCount, b_size);

  uintE* cts = newA(uintE, n_blocks+1);
  size_t* block_offs = newA(size_t, n_blocks+1);

  auto offsets_m = make_in_imap<uintT>(m, [&] (size_t i) { return offsets[i]; });
  auto lt = [] (const uintT& l, const uintT& r) { return l < r; };
  parallel_for(size_t i=0; i<n_blocks; i++) {
    size_t s_val = i*b_size;
    block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
  }
  block_offs[n_blocks] = m;
  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      // start and end are offsets in [m]
      size_t start = block_offs[i];
      size_t end = block_offs[i+1];
      uintT start_o = offsets[start];
      uintT k = start_o;
      for (size_t j=start; j<end; j++) {
        uintE v = indices.vtx(j);
        size_t num_in = frontierVertices[j].decodeOutNghSparseSeq(v, k, f, g);
        k += num_in;
      }
      cts[i] = (k - start_o);
    } else {
      cts[i] = 0;
    }
  }

  long outSize = sequence::plusScan(cts, cts, n_blocks);
  cts[n_blocks] = outSize;

  S* out = newA(S, outSize);

  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      size_t start = block_offs[i];
      size_t start_o = offsets[start];
      size_t out_off = cts[i];
      size_t block_size = cts[i+1] - out_off;
      for (size_t j=0; j<block_size; j++) {
        out[out_off + j] = outEdges[start_o + j];
      }
    }
  }
  free(outEdges); free(cts); free(block_offs);

  if (fl & remove_duplicates) {
    if (GA.flags == NULL) {
      GA.flags = newA(uintE, n);
      parallel_for(size_t i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
    }
    auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(out[i]); };
    remDuplicates(get_key, GA.flags, outSize, n);
    S* nextIndices = newA(S, outSize);
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
    free(out);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  }
  return vertexSubsetData<data>(n, outSize, out);
}
template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter(graph<vertex>& GA,
    vertex** frontierVertices, VS& indices, uintT* offsets, uintT m, F& f,
    const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  long outEdgeCount = sequence::plusScan(offsets, offsets, m);
  S* outEdges = newA(S, outEdgeCount);

  auto g = get_emsparse_no_filter_gen<data>(outEdges);

  // binary-search into scan to map workers->chunks
  size_t b_size = 10000;
  size_t n_blocks = nblocks(outEdgeCount, b_size);

  uintE* cts = newA(uintE, n_blocks+1);
  size_t* block_offs = newA(size_t, n_blocks+1);

  auto offsets_m = make_in_imap<uintT>(m, [&] (size_t i) { return offsets[i]; });
  auto lt = [] (const uintT& l, const uintT& r) { return l < r; };
  parallel_for(size_t i=0; i<n_blocks; i++) {
    size_t s_val = i*b_size;
    block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
  }
  block_offs[n_blocks] = m;
  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      // start and end are offsets in [m]
      size_t start = block_offs[i];
      size_t end = block_offs[i+1];
      uintT start_o = offsets[start];
      uintT k = start_o;
      for (size_t j=start; j<end; j++) {
        uintE v = indices.vtx(j);
        size_t num_in = frontierVertices[j]->decodeOutNghSparseSeq(v, k, f, g);
        k += num_in;
      }
      cts[i] = (k - start_o);
    } else {
      cts[i] = 0;
    }
  }

  long outSize = sequence::plusScan(cts, cts, n_blocks);
  cts[n_blocks] = outSize;

  S* out = newA(S, outSize);

  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      size_t start = block_offs[i];
      size_t start_o = offsets[start];
      size_t out_off = cts[i];
      size_t block_size = cts[i+1] - out_off;
      for (size_t j=0; j<block_size; j++) {
        out[out_off + j] = outEdges[start_o + j];
      }
    }
  }
  free(outEdges); free(cts); free(block_offs);

  if (fl & remove_duplicates) {
    if (GA.flags == NULL) {
      GA.flags = newA(uintE, n);
      parallel_for(size_t i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
    }
    auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(out[i]); };
    remDuplicates(get_key, GA.flags, outSize, n);
    S* nextIndices = newA(S, outSize);
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
    free(out);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  }
  return vertexSubsetData<data>(n, outSize, out);
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapData(graph<vertex>& GA, VS &vs, F f,
    intT threshold = -1, const flags& fl=0) {
  long numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  if(threshold == -1) threshold = numEdges/20; //default threshold
  // vertex *G = GA.getvertex();
  // ChronoTimer cter;
  if (numVertices != vs.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  if (vs.size() == 0) return vertexSubsetData<data>(numVertices);
  vs.toSparse();
  uintT* degrees = newA(uintT, m);
  /*
  vertex* frontierVertices = newA(vertex,m);
  // vertex * frontierVertices = new vertex[m];
  {parallel_for (size_t i=0; i < m; i++) {
    uintE v_id = vs.vtx(i);
    vertex v = G[v_id];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
  }}
  */

  vertex ** frontierVertices = newA(vertex*, m);
  {parallel_for (size_t i=0; i<m; i++) {
    uintE v_id = vs.vtx(i);
    degrees[i] = GA.getvertex(v_id)->getOutDegree();
    frontierVertices[i] = GA.getvertex(v_id);
  }}

  uintT outDegrees = sequence::plusReduce(degrees, m);
  // cout << cter.elapsed_milli() << endl;
  // cout << m + outDegrees << " > " << threshold << (fl&dense_forward) << endl;
  if (outDegrees == 0) return vertexSubsetData<data>(numVertices);
  if (m + outDegrees > threshold) {
    vs.toDense();
    free(degrees);
    // free(frontierVertices);
    delete[] frontierVertices;
    return (fl & dense_forward) ?
      edgeMapDenseForward<data, vertex, VS, F>(&GA, vs, f, fl) :
      edgeMapDense<data, vertex, VS, F>(&GA, vs, f, fl);
  } else {
    auto vs_out =
      (should_output(fl) && fl & sparse_no_filter) ? // only call snof when we output
      edgeMapSparse_no_filter<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl) :
      edgeMapSparse<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl);
    free(degrees); free(frontierVertices);
    return vs_out;
  }
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class vertex, class VS, class F>
vertexSubset edgeMap(graph<vertex> * GA, VS& vs, F f,
    intT threshold = -1, const flags& fl=0) {
  return edgeMapData<pbbs::empty>(*GA, vs, f, threshold, fl);
}

// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
// Weighted graphs are not yet supported, but this should be easy to do.
template <class vertex, class P>
vertexSubsetData<uintE> packEdges(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  using S = tuple<uintE, uintE>;
  vs.toSparse();
  // vertex* G = GA.getvertex(); 
  long m = vs.numNonzeros(); long n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto degrees = array_imap<uintT>(m);
  granular_for(i, 0, m, (m > 2000), {
    uintE v = vs.vtx(i);
    degrees[i] = GA.getvertex(v)->getOutDegree();
  });
  long outEdgeCount = pbbs::scan_add(degrees, degrees);
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }

  bool* bits = newA(bool, outEdgeCount);
  uintE* tmp1 = newA(uintE, outEdgeCount);
  uintE* tmp2 = newA(uintE, outEdgeCount);
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = GA.getvertex(v)->packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = GA.getvertex(v)->packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
    }
  }
  free(bits); free(tmp1); free(tmp2);
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}

template <class vertex, class P>
vertexSubsetData<uintE> edgeMapFilter(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  vs.toSparse();
  if (fl & pack_edges) {
    return packEdges<vertex, P>(GA, vs, p, fl);
  }
  // vertex* G = GA.getvertex(); 
  long m = vs.numNonzeros(); long n = vs.numRows();
  using S = tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = GA.getvertex(v)->countOutNgh(v, p);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = GA.getvertex(v)->countOutNgh(v, p);
    }
  }
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}



//*****VERTEX FUNCTIONS*****

template <class F, class VS, typename std::enable_if<
  !std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i), V.vtxData(i));
    }
  }
}

template <class VS, class F, typename std::enable_if<
  std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i);
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i));
    }
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

template <class F>
vertexSubset vertexFilter2(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  {parallel_for(size_t i=0; i<m; i++) {
    uintE v = V.vtx(i);
    bits[i] = filter(v);
  }}
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}

template <class data, class F>
vertexSubset vertexFilter2(vertexSubsetData<data> V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  parallel_for(size_t i=0; i<m; i++) {
    auto t = V.vtxAndData(i);
    bits[i] = filter(std::get<0>(t), std::get<1>(t));
  }
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}


//cond function that always returns true
inline bool cond_true (intT d) { return 1; }

template<class vertex>
void Compute(graph<vertex>&, commandLine);

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  bool hybrid_vertex = P.getOptionValue("-h");
  bool version_graph = P.getOptionValue("-g");
  int delta_num = P.getOptionIntValue("-n", 9);
  int delta_type = P.getOptionIntValue("-type", 0);
  double add_rate = P.getOptionDoubleValue("-r", 2);
  double delta_rate = P.getOptionDoubleValue("-f", 2);
  int offset_wall = P.getOptionIntValue("-w", 0);

  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",20);

  string exe_file = split_last(argv[0], '/');
  int pid = get_pid(exe_file.c_str());
  cout << exe_file << " " << pid << endl;

  cout.setf(ios::fixed, ios::floatfield);

  string filename(iFile);
  // string basedir = get_basedir(filename, delta_rate, add_rate, delta_type);
  string basedir = "/public/home/tangwei/graphdata/flickr/finaltreetest/";

  /////  code for generate delta
{

  /* 
  int delta_step = 1;
  graph<symmetricVertex> G = readGraphFromFile<symmetricVertex>(iFile);
  vector<delta_log<symmetricVertex>> deltalogs;
  vector<delta<symmetricVertex>> deltas;
  // graph<symmetricVertex> ga = readGraphFromFile<symmetricVertex>(iFile);
  string savedir = "/public/home/tangwei/graphdata/flickr/finaltree/";
  // string savedir;
  /* 
  if (delta_rate == 2) {
    delta_rate = 0.001;
  }
  if (add_rate == 2) {
    add_rate = 0.9;
  }
  if (filename.find("flickr") != filename.npos) {
    ostringstream str1, str2;
    str1 << delta_rate;
    string delta_rate_string = str1.str();
    str2 << add_rate;
    string add_rate_string = str2.str();
    savedir = string("/public/home/tangwei/graphdata/flickr/d") + delta_rate_string + "-a" + add_rate_string + '/';

    // savedir = "/public/home/tangwei/graphdata/flickr/d" + to_string(delta_rate) + "-a" + to_string(add_rate) + '/'; 
    // savedir = "/public/home/tangwei/graphdata/flickr/deletedelta/";
  } else if (filename.find("enwiki") != filename.npos) {
    savedir = "/public/home/tangwei/graphdata/enwiki/smalldelta/";
  } else if (filename.find("orkut") != filename.npos) {
    savedir = "/public/home/tangwei/graphdata/orkut/smalldelta/";
  } else if (filename.find("journal") != filename.npos) {
    savedir = "/public/home/tangwei/graphdata/livejournal/smalldelta/";
  } else if (filename.find("dewiki") != filename.npos) {
    savedir = "/public/home/tangwei/graphdata/link-dynamic-dewiki/smalldelta/";
  } else if (filename.find("twitter") != filename.npos) {
    savedir = "/public/home/tangwei/graphdata/twitter/smalldelta/";
  } else {
    cout << "file " << iFile << " dont have related delta file, please reach admin for help" << endl;
    return 0;
  }
  *
  // 
  
  for (auto i=1; i<=delta_num; i++) {
    uintT new_from = get_new_from(i-1, delta_step);
    cout << new_from << " to " << i << endl;
    vector<uintT> path = get_revert_path(deltas, G.get_version(), new_from);
    for (auto j : path) {
      auto tmp = deltas[j];
      if (tmp.ve == G.get_version()) {
        revert(G, tmp);
      } else if (tmp.vs == G.get_version()) {
        apply(G, tmp);
      } else {
        cout << "error! delta in revert path not valid" << endl;
      }
    }
    cout << "finish get path" << endl;
    delta_log<symmetricVertex> dlg = delta_log<symmetricVertex>(G, 0.9, 0.001, true);
    cout << "finish get dlg" << endl;
    cout << dlg.delta_rate << endl;
    cout << dlg.deltaLog.size() << endl;
    dlg.ver_end = i;
    delta<symmetricVertex> dlt = delta<symmetricVertex>(dlg, G);
    cout << "finish get dlt" << endl;
    deltalogs.push_back(dlg);
    deltas.push_back(dlt);
    // apply(G, dlt);
    cout << "dirty after roll " << G.accessAllEdges() << endl;

    string savefile = savedir + to_string(i-1) + '-' + to_string(i);
    if (delta_step == 0) {
      dlt.write_edge_entry(savefile.c_str());
      cout << "finish write" << endl;
    }
    uintT s = i-1;
    uintT e = i;
    path = get_revert_path(deltas, s, e);
    vector<delta<symmetricVertex> *> ppath;
    for (auto & p: path) {
      ppath.push_back(&deltas[p]);
      cout << deltas[p].vs << " " << deltas[p].ve << endl;
    }
    delta<symmetricVertex> merged = delta<symmetricVertex>(G, ppath, s, e);
    merged.vs = s;
    merged.ve = e;
    
    merged.write_edge_entry(savefile.c_str());
    cout << "finish write " << s << endl;
    apply(G, merged);
  }
  G.del();
  // ga.del();
  return 0;
  // */
}
  /////  end of code for generate tree delta
  ChronoTimer cter;
  if (!hybrid_vertex) {
    startTime();
    graph<HVertex> G = readHGraphFromFile(iFile, false);
    // graph<symmetricVertex> G = readGraphFromFile<symmetricVertex>(iFile);
      // readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    nextTime("load graph ");
    print_mem(pid);

    if (!version_graph) {
      // bigDelta<symmetricVertex> bdelta;
      bigDelta<HVertex> bdelta;

      // get delta and merge into a bigDelta
      for (auto i = 0; i < delta_num; i++) {
        cout << "delta : " << i << "-" << i+1 << endl;
        string filename = basedir + to_string(i) + '-' + to_string(i+1);
        // delta_log<symmetricVertex> dlg = load_deltalog_from_file<symmetricVertex>(G, filename.c_str());
        // delta<symmetricVertex> dlt = delta<symmetricVertex>(dlg, G);
        delta_log<HVertex> dlg = load_deltalog_from_file<HVertex>(G, filename.c_str());
        delta<HVertex> dlt = delta<HVertex>(dlg, G);
        bdelta.append(dlt);
        forward(G, bdelta);
        if (!(i%5)) {
          print_mem(pid);
          // startTime();
          cter.start();
          cout << "dirty data " << G.accessAllEdges() << endl;
          // nextTime("access all edge");
          cout << "access all edge : " << cter.elapsed() << endl;
        }
      }

      // switch and compute every 5 version
      for(int r = 0; r < delta_num; r++) {
        if (r%5) { continue; }
        jump(G, bdelta, delta_num);
        string print_data = "switch time to " + to_string(r) + " : ";
        // startTime();
        cter.start();
        jump(G, bdelta, r);
        // nextTime("access all edge");
        cout << "switch time to " + to_string(r) + " " << cter.elapsed() << endl;
      }
      return 0;
      
      // random switch with skew in version
      for(int r=0; r<distOrder::get_size(); r++) {
        int randi = distOrder::get(r, delta_num);
        string print_data1 = "addtime1 to version " + to_string(r) + " : ";
        string print_data2 = "addtime2 to version " + to_string(r) + " : ";
        jump(G, bdelta, delta_num);
        startTime();
        jump(G, bdelta, randi);
        nextTime(print_data1)
        Compute(G,P);
        nextTime(print_data2);
      }
    } else {
      // versionGraph<symmetricVertex> vg;
      versionGraph<HVertex> vg;

      // get delta and merge into a version graph
      for (auto i=0; i < delta_num; i++) {
        cout << "version graph delta : " << i << "-" << i+1 << endl;
        string filename = basedir + to_string(i) + '-' + to_string(i+1);
        // delta_log<symmetricVertex> dlg = load_deltalog_from_file<symmetricVertex>(G, filename.c_str());
        // delta<symmetricVertex> *dlt = new delta<symmetricVertex>(dlg, G);
        delta_log<HVertex> dlg = load_deltalog_from_file<HVertex>(G, filename.c_str());
        delta<HVertex> * dlt = new delta<HVertex>(dlg, G);
        uintT tmp = dlt->ve;

        vg.append(G, dlt);

        jump(G, vg, tmp);
        vg.shrink();
        G.shrink();
        if (!(i%5)) {
          print_mem(pid);
          // startTime();
          cter.start();
          cout << "dirty data " << G.accessAllEdges() << endl;
          // nextTime("access all edge");
          cout << "access all edge : " << cter.elapsed() << endl;
        }
      }

      // switch and compute every 5 version
      for(int r = 0; r < delta_num; r++) {
        if (r%5) { continue; }
        jump(G, vg, delta_num);
        string print_data = "switch time to " + to_string(r) + " : ";
        // startTime();
        cter.start();
        jump(G, vg, r);
        // nextTime("access all edge");
        cout << "switch time to " + to_string(r) + " " << cter.elapsed() << endl;
      }
      return 0;

      // random switch with skew in version
      for(int r=0; r<distOrder::get_size(); r++) {
        int randi = distOrder::get(r, delta_num);
        string print_data1 = "addtime1 to version " + to_string(r) + " : ";
        string print_data2 = "addtime2 to version " + to_string(r) + " : ";
        jump(G, vg, delta_num);
        startTime();
        jump(G, vg, randi);
        nextTime(print_data1);
        Compute(G,P);
        nextTime(print_data2);
      }
    }
    G.del();
  } else { // code for hybrid vertex
    cout << "vertex type : hybrid " << endl;
    startTime();
    graph<HVertex> G = readHGraphFromFile(iFile, true, offset_wall);
    nextTime("load graph ");
    // startTime();

    if (!version_graph) {
      bigDelta<HVertex> bdelta;

      // get delta and merge into a bigDelta
      for (auto i = 0; i < delta_num; i++) {
        cout << "delta : " << i << "-" << i+1 << endl;
        string filename = basedir + to_string(i) + '-' + to_string(i+1);
        delta_log<HVertex> dlg = load_deltalog_from_file<HVertex>(G, filename.c_str());
        delta<HVertex> dlt = delta<HVertex>(dlg, G);
        bdelta.append(dlt);
        forward(G, bdelta);
        if (!(i%5)) {
          print_mem(pid);
          // startTime();
          cter.start();
          cout << "dirty data " << G.accessAllEdges() << endl;
          // nextTime("access all edge");
          cout << "access all edge : " << cter.elapsed() << endl;
        }
      }

      // switch and compute every 5 version
      for(int r = 0; r < delta_num; r++) {
        if (r%5) { continue; }
        jump(G, bdelta, delta_num);
        string print_data = "switch time to " + to_string(r) + " : ";
        // startTime();
        cter.start();
        jump(G, bdelta, r);
        // nextTime("access all edge");
        cout << "switch time to " + to_string(r) + " " << cter.elapsed() << endl;
      }

      return 0;

      // random switch with skew in version
      for(int r=0; r<distOrder::get_size(); r++) {
        int randi = distOrder::get(r, delta_num);
        string print_data1 = "addtime1 to version " + to_string(r) + " : ";
        string print_data2 = "addtime2 to version " + to_string(r) + " : ";
        jump(G, bdelta, delta_num);
        startTime();
        jump(G, bdelta, randi);
        nextTime(print_data1);
        Compute(G,P);
        nextTime(print_data2);
      }
    } else {
      versionGraph<HVertex> vg;
      // get delta and merge into a version graph
      for (auto i=0; i<delta_num; i++) {
        cout << "version graph delta : " << i << "-" << i+1 << endl;
        string filename = basedir + to_string(i) + '-' + to_string(i+1);
        // print_mem(pid);
        delta_log<HVertex> dlg = load_deltalog_from_file<HVertex>(G, filename.c_str());
        // print_mem(pid);
        delta<HVertex> *dlt = new delta<HVertex>(dlg, G);
        // print_mem(pid);
        uintT tmp = dlt->ve;
        vg.append(G, dlt);
        // print_mem(pid);
        vg.shrink();
        G.shrink();
        jump(G, vg, tmp);
        if (!(i%5)) {
          print_mem(pid);
          // startTime();
          cter.start();
          cout << "dirty data " << G.accessAllEdges() << endl;
          // nextTime("access all edge");
          cout << "access all edge : " << cter.elapsed() << endl;
        }
      }

      // switch and compute every 5 version
      for(int r = 0; r < delta_num; r++) {
        if (r%5) { continue; }
        jump(G, vg, delta_num);
        string print_data = "switch time to " + to_string(r) + " : ";
        // startTime();
        cter.start();
        jump(G, vg, r);
        // nextTime("access all edge");
        cout << "switch time to " + to_string(r) + " " << cter.elapsed() << endl;
      }
      return 0;

      // random switch with skew in version
      for(int r=0; r<distOrder::get_size(); r++) {
        int randi = distOrder::get(r, delta_num);
        string print_data1 = "addtime1 to version " + to_string(r) + " : ";
        string print_data2 = "addtime2 to version " + to_string(r) + " : ";
        jump(G, vg, delta_num);
        startTime();
        jump(G, vg, randi);
        nextTime(print_data1);
        Compute(G,P);
        nextTime(print_data2);
      }
    }
    G.del();
  }
  return 0;
}

#endif
