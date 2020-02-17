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
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include "parallel.h"
#include "blockRadixSort.h"
#include "quickSort.h"
#include "utils.h"
#include "graph.h"
// #include "vector.h"
#include <vector>
#include "get_mem.h"
#include "delta.h"
#include "myutil.h"
using namespace std;

typedef pair<uintE,uintE> intPair;
typedef pair<uintE, pair<uintE,intE> > intTriple;

template <class E>
struct pairFirstCmp {
  bool operator() (pair<uintE,E> a, pair<uintE,E> b) {
    return a.first < b.first; }
};

template <class E>
struct getFirst {uintE operator() (pair<uintE,E> a) {return a.first;} };

template <class IntType>
struct pairBothCmp {
  bool operator() (pair<uintE,IntType> a, pair<uintE,IntType> b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  }
};

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long n; // total number of characters
  char* Chars;  // array storing all strings
  long m; // number of substrings
  char** Strings; // pointers to strings (all should be null terminated)
  words() {}
words(char* C, long nn, char** S, long mm)
: Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {free(Chars); free(Strings);}
};

inline bool isSpace(char c) {
  switch (c)  {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ' : return true;
  default : return false;
  }
}

_seq<char> mmapStringFromFile(const char *filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG (sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char *p = static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
//  char *bytes = newA(char, n);
//  parallel_for(size_t i=0; i<n; i++) {
//    bytes[i] = p[i];
//  }
//  if (munmap(p, sb.st_size) == -1) {
//    perror("munmap");
//    exit(-1);
//  }
//  cout << "mmapped" << endl;
//  free(bytes);
//  exit(0);
  return _seq<char>(p, n);
}

_seq<char> readStringFromFile(char *fileName) {
  ifstream file (fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg (0, ios::beg);
  long n = end - file.tellg();
  char* bytes = newA(char,n+1);
  file.read (bytes,n);
  file.close();
  return _seq<char>(bytes,n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, long n) {
  {parallel_for (long i=0; i < n; i++)
      if (isSpace(Str[i])) Str[i] = 0; }

  // mark start of words
  bool *FL = newA(bool,n);
  FL[0] = Str[0];
  {parallel_for (long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];}

  // offset for each start of word
  _seq<long> Off = sequence::packIndex<long>(FL, n);
  long m = Off.n;
  long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char*, m);
  {parallel_for (long j=0; j < m; j++) SA[j] = Str+offsets[j];}

  free(offsets); free(FL);
  return words(Str,n,SA,m);
}

template <class vertex>
delta<vertex> readDeltaFromFile(char* fname) {
  words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);

  if (W.Strings[0] != (string) "AdjacencyGraphDelta") {
    cout << "bad input file" << endl;
    abort();
  }

  int version_start = atoi(W.Strings[1]);
  int version_end = atoi(W.Strings[2]);

  long len = W.m - 1;
  long len_v = atol(W.Strings[3]);
  long len_p = atol(W.Strings[4]);
  long len_d = atol(W.Strings[5]);
  if (len != len_v + len_p + len_d) {
    cout << "bad input file" << endl;
    abort();
  }

  vector<uintE> vertex_delta;
  vector<uintE> pos_delta;
  vector<uintE> dst_delta;
  vertex_delta.reserve(len_v);
  pos_delta.reserve(len_p);
  dst_delta.reserve(len_d);
  {parallel_for(long i=0; i < len_v; i++) vertex_delta[i, atol(W.Strings[i+6])];}
  {parallel_for(long j=0; j < len_p; j++) pos_delta[j, atol(W.Strings[j+len_v+6])];}
  {parallel_for(long k=0; k < len_d; k++) dst_delta[k, atol(W.Strings[k+len_v+len_p+6])];}

  return delta<vertex>(version_start, version_end, vertex_delta, pos_delta, dst_delta);
}

template <class vertex>
delta<vertex> readDeltaFromLog(string fname, graph<vertex> & ga) {
  ifstream fin(fname);
  if (!fin) {
    cout << "bad input delta file " << endl;
  }

  cout << "begin to parse delta file" << endl;

  string line;
  vector<intTriple> logs;
  int srcv, dstv, wgh;
  while(getline(fin, line)) {
    if (line[0] == '#' || line[0] == '%') {
      continue;
    }
    istringstream iss(line);
    iss >> srcv >> dstv >> wgh;
    logs.push_back(make_pair(srcv, make_pair(dstv, wgh)));
  }
  return delta<vertex>(ga, logs);
}

template <class vertex>
graph<vertex> readGraphFromFile(char *fname) {
  words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);
  if (W.Strings[0] != (string) "AdjacencyGraph") {
    cout << "Bad input file : bad head" << endl;
    cout << W.Strings[0] << endl;
    exit(0);
  }

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);

  if (len != n + m + 2) {
    cout << "Bad input file : bad size in head" << endl;
    cout << "len: " << len << " n: " << n << " m: " << m << endl;
    cout << W.Strings[0] << W.Strings[len] << endl;
    exit(0);
  }

  uintT* offsets = newA(uintT,n);
  uintE* edges = newA(uintE,m);
  {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long i=0; i<m; i++) {
    edges[i] = atol(W.Strings[i+n+3]);
  }}

  W.del();

  cout << "mid of read graph with n " << n << " and m " << m << endl;

  vector<vertex> v;
  cout << v.capacity() << endl;
  vertex vtmp;
  v.resize(n, vtmp);
  cout << v.capacity() << endl;

  {parallel_for(uintT i = 0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = (i == n - 1) ? m : offsets[i + 1];
    v[i].reserve(l-o);
    for (uintT j = o; j < l; j++)
      v[i].push_back(edges[j]);
  }}//Put the outneighbors of each edge into the corresponding vector


  free(offsets);
  free(edges);

  // Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,NULL);
  Empty_Mem * mem = new Empty_Mem();
  return graph<vertex>(v,n,m,mem);
}

uintT get_offset_wall(uintT * arr, long number, long length, float percent=0.05) {
  vector<uintT> v;
  
  for (auto i=0; i<number; i++) {
    if (i==number-1) {
      v.push_back(length-arr[i]);
    } else {
      v.push_back(arr[i+1] - arr[i]);
    }
  }
  
  sort(v.begin(), v.end());
  if (percent <=0) {
    return v[number-1];
  } else if (percent >= 1) {
    return v[0];
  }
  return v[int(number*(1-percent))];
}

graph<HVertex> readHGraphFromFile(char *fname, bool hybrid=true, int _offset_wall=0) {
    words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);
  if (W.Strings[0] != (string) "AdjacencyGraph") {
    cout << "Bad input file : bad head" << endl;
    cout << W.Strings[0] << endl;
    exit(0);
  }

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);

  if (len != n + m + 2) {
    cout << "Bad input file : bad size in head" << endl;
    cout << "len: " << len << " n: " << n << " m: " << m << endl;
    cout << W.Strings[0] << W.Strings[len] << endl;
    exit(0);
  }

  uintT* offsets = newA(uintT,n);
  uintE* edges = newA(uintE,m);
  {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long i=0; i<m; i++) {
    edges[i] = atol(W.Strings[i+n+3]);
  }}
  
  W.del();
  
  // uintT offset_wall = hybrid ? get_offset_wall(offsets, n, m) : 0;
  uintT offset_wall;
  if (hybrid) {
    if (_offset_wall == 0) {
      offset_wall = get_offset_wall(offsets, n, m);
    } else {
      offset_wall = static_cast<uintT>(_offset_wall);
    }
  } else {
    offset_wall = 0;
  }

  cout << "get wall " << offset_wall << endl;
  vector<HVertex> v;

  // unordered_map<uintE, bool> hdv;
  v.reserve(n);
  {parallel_for(uintT i=0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = (i == n - 1) ? m : offsets[i + 1];
    v.emplace_back(edges+o, l-o, !hybrid ||  (l-o > offset_wall));
  }}

  free(offsets);
  free(edges);
  Empty_Mem * mem = new Empty_Mem();

  // cout << "before return" << endl;
  return graph<HVertex>(v,n,m,mem, hybrid);
}

template <class vertex>
graph<vertex> readGraphFromFile(char* fname, bool isSymmetric, bool mmap) {
  words W;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    // Cannot mutate the graph unless we copy.
    parallel_for(size_t i=0; i<S.n; i++) {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    _seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n); // you dont free S.A here. I did and I failed
  }
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyGraph") {
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph") {
#endif
    cout << "Bad input file : bad head" << endl;
    cout << W.Strings[0] << endl;
    exit(0);
  }

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
#ifndef WEIGHTED
  if (len != n + m + 2) {
#else
  if (len != n + 2*m + 2) {
#endif
    cout << "Bad input file : bad size in head" << endl;
    cout << "len: " << len << " n: " << n << " m: " << m << endl;
    cout << W.Strings[0] << W.Strings[len] << endl;
    exit(0);
  }

  uintT* offsets = newA(uintT,n);
#ifndef WEIGHTED
  uintE* edges = newA(uintE,m);
#else
  intE* edges = newA(intE,2*m);
#endif

  {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long i=0; i<m; i++) {
#ifndef WEIGHTED
      edges[i] = atol(W.Strings[i+n+3]);
#else
      edges[2*i] = atol(W.Strings[i+n+3]);
      edges[2*i+1] = atol(W.Strings[i+n+m+3]);
#endif
    }}
  W.del(); // to deal with performance bug in malloc

  vector<vertex> v ;
  v.resize(n);
//
#ifndef WEIGHTED
  {parallel_for(uintT i = 0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = (i == n - 1) ? m : offsets[i + 1];
    v[i].outNeighbors.reserve(l-o);
    for (uintT j = o; j < l; j++)
      v[i].outNeighbors.push_back(edges[j]);
  }}//Put the outneighbors of each edge into the corresponding vector
#endif

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    {parallel_for(long i=0; i<n; i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair, m);
#else
    intTriple* temp = newA(intTriple, m);
#endif
    {parallel_for(long i=0;i<n;i++){
      uintT o = offsets[i];
      for(uintT j=0;j < v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsets);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<uintE>());//sort temp by the first
#else
    quickSort(temp,m,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<intPair>());
#else
    quickSort(temp,m,pairFirstCmp<intPair>());
#endif
#endif

    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE,2*m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {parallel_for(long i=1;i<m;i++) {
#ifndef WEIGHTED
      inEdges[i] = temp[i].second;
#else
      inEdges[2*i] = temp[i].second.first;
      inEdges[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);
#ifndef WEIGHTED
	{parallel_for(uintT i = 0; i < n; i++) {
		uintT o = tOffsets[i];
		uintT l = (i == n - 1) ? m : tOffsets[i + 1];
    v[i].outNeighbors.reserve(l-o);
		for (uintT j = o; j < l; j++)
			v[i].outNeighbors.push_back(inEdges[j]);
	}}//Put the inneighbors of each edge into the corresponding vector
#else
#endif
    

    free(tOffsets);
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges,inEdges);
    return graph<vertex>(v,n,m,mem);
  }
  else {
    free(offsets);
    free(edges);
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,NULL);
    return graph<vertex>(v,n,m,mem);
  }
}

template <class vertex>
graph<vertex> readGraphFromBinary(char* iFile, bool isSymmetric) {
  char* config = (char*) ".config";
  char* adj = (char*) ".adj";
  char* idx = (char*) ".idx";
  char configFile[strlen(iFile)+strlen(config)+1];
  char adjFile[strlen(iFile)+strlen(adj)+1];
  char idxFile[strlen(iFile)+strlen(idx)+1];
  *configFile = *adjFile = *idxFile = '\0';
  strcat(configFile,iFile);
  strcat(adjFile,iFile);
  strcat(idxFile,iFile);
  strcat(configFile,config);
  strcat(adjFile,adj);
  strcat(idxFile,idx);

  ifstream in(configFile, ifstream::in);
  long n;
  in >> n;
  in.close();

  ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints
  in2.seekg(0, ios::end);
  long size = in2.tellg();
  in2.seekg(0);
#ifdef WEIGHTED
  long m = size/(2*sizeof(uint));
#else
  long m = size/sizeof(uint);
#endif
  char* s = (char *) malloc(size);
  in2.read(s,size);
  in2.close();
  uintE* edges = (uintE*) s;

  ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if(n != size/sizeof(intT)) { cout << "File size wrong\n"; abort(); }

  char* t = (char *) malloc(size);
  in3.read(t,size);
  in3.close();
  uintT* offsets = (uintT*) t;

  vertex* v = newA(vertex,n);
#ifdef WEIGHTED
  intE* edgesAndWeights = newA(intE,2*m);
  {parallel_for(long i=0;i<m;i++) {
    edgesAndWeights[2*i] = edges[i];
    edgesAndWeights[2*i+1] = edges[i+m];
    }}
  //free(edges);
#endif
  {parallel_for(long i=0;i<n;i++) {
    uintT o = offsets[i];
    uintT l = ((i==n-1) ? m : offsets[i+1])-offsets[i];
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      v[i].setOutNeighbors((uintE*)edges+o);
#else
      v[i].setOutNeighbors(edgesAndWeights+2*o);
#endif
    }}

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair,m);
#else
    intTriple* temp = newA(intTriple,m);
#endif
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsets);
#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<uintE>());
#else
    quickSort(temp,m,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<intPair>());
#else
    quickSort(temp,m,pairFirstCmp<intPair>());
#endif
#endif
    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE,2*m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {parallel_for(long i=1;i<m;i++) {
#ifndef WEIGHTED
      inEdges[i] = temp[i].second;
#else
      inEdges[2*i] = temp[i].second.first;
      inEdges[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    free(temp);
    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);
    {parallel_for(long i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
#ifndef WEIGHTED
      v[i].setInNeighbors((uintE*)inEdges+o);
#else
      v[i].setInNeighbors((intE*)(inEdges+2*o));
#endif
      }}
    free(tOffsets);
#ifndef WEIGHTED
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges,inEdges);
    return graph<vertex>(v,n,m,mem);
#else
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights,inEdges);
    return graph<vertex>(v,n,m,mem);
#endif
  }
  free(offsets);
#ifndef WEIGHTED
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges);
  return graph<vertex>(v,n,m,mem);
#else
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights);
  return graph<vertex>(v,n,m,mem);
#endif
}

template <class vertex>
graph<vertex> readGraph(char* iFile, bool compressed, bool symmetric, bool binary, bool mmap) {
  if(binary) return readGraphFromBinary<vertex>(iFile,symmetric);
  // else return readGraphFromFile<vertex>(iFile,symmetric,mmap);
  else return readGraphFromFile<vertex>(iFile);
}

// template <class vertex>
// graph<vertex> readCompressedGraph(char* fname, bool isSymmetric, bool mmap) {
//   char* s;
//   if (mmap) {
//     _seq<char> S = mmapStringFromFile(fname);
//     // Cannot mutate graph unless we copy.
//     char *bytes = newA(char, S.n);
//     parallel_for(size_t i=0; i<S.n; i++) {
//       bytes[i] = S.A[i];
//     }
//     if (munmap(S.A, S.n) == -1) {
//       perror("munmap");
//       exit(-1);
//     }
//     s = bytes;
//   } else {
//     ifstream in(fname,ifstream::in |ios::binary);
//     in.seekg(0,ios::end);
//     long size = in.tellg();
//     in.seekg(0);
//     cout << "size = " << size << endl;
//     s = (char*) malloc(size);
//     in.read(s,size);
//     in.close();
//   }

//   long* sizes = (long*) s;
//   long n = sizes[0], m = sizes[1], totalSpace = sizes[2];

//   cout << "n = "<<n<<" m = "<<m<<" totalSpace = "<<totalSpace<<endl;
//   cout << "reading file..."<<endl;

//   uintT* offsets = (uintT*) (s+3*sizeof(long));
//   long skip = 3*sizeof(long) + (n+1)*sizeof(intT);
//   uintE* Degrees = (uintE*) (s+skip);
//   skip+= n*sizeof(intE);
//   uchar* edges = (uchar*)(s+skip);

//   uintT* inOffsets;
//   uchar* inEdges;
//   uintE* inDegrees;
//   if(!isSymmetric){
//     skip += totalSpace;
//     uchar* inData = (uchar*)(s + skip);
//     sizes = (long*) inData;
//     long inTotalSpace = sizes[0];
//     cout << "inTotalSpace = "<<inTotalSpace<<endl;
//     skip += sizeof(long);
//     inOffsets = (uintT*) (s + skip);
//     skip += (n+1)*sizeof(uintT);
//     inDegrees = (uintE*)(s+skip);
//     skip += n*sizeof(uintE);
//     inEdges = (uchar*)(s + skip);
//   } else {
//     inOffsets = offsets;
//     inEdges = edges;
//     inDegrees = Degrees;
//   }


//   vertex *V = newA(vertex,n);
//   parallel_for(long i=0;i<n;i++) {
//     long o = offsets[i];
//     uintT d = Degrees[i];
//     V[i].setOutDegree(d);
//     V[i].setOutNeighbors(edges+o);
//   }

//   if(sizeof(vertex) == sizeof(compressedAsymmetricVertex)){
//     parallel_for(long i=0;i<n;i++) {
//       long o = inOffsets[i];
//       uintT d = inDegrees[i];
//       V[i].setInDegree(d);
//       V[i].setInNeighbors(inEdges+o);
//     }
//   }

//   cout << "creating graph..."<<endl;
//   Compressed_Mem<vertex>* mem = new Compressed_Mem<vertex>(V, s);

//   graph<vertex> G(V,n,m,mem);
//   return G;
// }
