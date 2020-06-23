#ifndef DELTA_H
#define DELTA_H

#include "graph.h"
#include <vector>
#include "vertex.h"
#include "quickSort.h"
#include <cstdlib>
#include <get_mem.h>
#include <unordered_map>
#include "delta_log.h"


bool equal(pair<uintE, uintE> a, pair<uintE, uintE> b);
bool lessthan(pair<uintE, uintE> a, pair<uintE, uintE> b);

// Pensieve delta, csr-like delta, with pos of delete entry.
template <class vertex>
struct delta {
  int vs, ve; // version start and version end
  vector<uintE> vertexs;
  vector<uintE> positions;
  vector<uintE> dstAndPos;
  int del_size = -1;
  delta() {}
  delta(int _vs, int _ve, vector<uintE>& v, vector<uintE>& p, vector<uintE>& dap) :
    vs(_vs), ve(_ve), vertexs(v), positions(p), dstAndPos(dap) {}


  delta(delta_log<vertex> &dlg, graph<vertex> & graph, bool origin = true) {
    if (graph.get_version() != dlg.ver) {
      cout  << "error in generate delta from deltalog, with version " 
            << dlg.ver << " -> " << dlg.ver_end 
            << " graph version " << graph.get_version() << endl; 
      abort();
    }

    gen_hybrid(dlg, graph, origin);
    vs = graph.get_version();
    ve = dlg.ver_end;
    del_size = -1;
    del_size = get_del_number();
  }

  delta(graph<vertex> & graph, vector<delta<vertex> *> path, uintT ver_begin, uintT ver_end) {
    vector<pair<uintE, uintE>> tmp_add, tmp_del;
    // cout << "into delta gen" << endl;
    for (auto i : path) {
      if (i->vs == ver_begin) {
        i->pack_add_items(&tmp_add);
        i->pack_del_items(&tmp_del);
        ver_begin = i->ve;
      } else if (i->ve == ver_begin) {
        i->pack_add_items(&tmp_del);
        i->pack_del_items(&tmp_add);
        ver_begin = i->vs;
      } else {
        cout << "error in trace path" << endl;
      }
    }

    if (ver_begin != ver_end) {
      cout << "may be delta in path not merged" << endl;
    }

    // 
    sort(tmp_add.begin(), tmp_add.end(), logLT());
    sort(tmp_del.begin(), tmp_del.end(), logLT());

    // drop items in both add and del.
    // save others to ret_add and ret_del.
    vector<pair<uintE, uintE>> ret_add, ret_del;
    uintT pa=0, pd=0;
    // using double pointer
    while (pa < tmp_add.size() && pd < tmp_del.size()) {
      if (equal(tmp_add[pa], tmp_del[pd])) {
        pa ++; pd++;
      } else if (lessthan(tmp_add[pa], tmp_del[pd])) {
        ret_add.push_back(make_pair(tmp_add[pa].first, tmp_add[pa].second));
        pa ++;
      } else {
        ret_del.push_back(make_pair(tmp_del[pd].first, tmp_del[pd].second));
        pd++;
      }
    }
    // process terminate case.
    if (pa == tmp_add.size()) {
      ret_del.insert(ret_del.end(), tmp_del.begin() + pd, tmp_del.end());
    } else if (pd == tmp_del.size()) {
      ret_add.insert(ret_add.end(), tmp_add.begin() + pa, tmp_add.end());
    }

    // cout << "while merge : add " << ret_add.size() << " and del " << ret_del.size() << endl;

    // use ret_add and ret_del to gen a delta_log on graph.
    delta_log<vertex> dlg = delta_log<vertex>(graph, ret_add, ret_del);

    // cout << "before gen " << endl;
    gen_hybrid(dlg, graph, false);

    vs = graph.get_version();
    // ve = get_path_end(path, ver_end);
    ve = ver_end;
    del_size = -1;
    del_size = get_del_number();
  }

  delta(const char * filename) {
    ifstream infile;
    infile.open(filename);
    if (!infile.is_open()) {
      cout << "fail to open file when try to get delta." << endl;
      cout << filename << endl;
      abort();
    }
    string header;
    infile >> header;
    if (header != "PENSIEVE_DELTA_FILE") {
      cout << "bad format. expect pensieve delta file, find " << header << endl;
    }
    infile >> vs >> ve;
    int a, b, c;
    infile >> a >> b >> c;
    vertexs.resize(a);
    positions.resize(b);
    dstAndPos.resize(c);
    for (auto i=0; i<a; i++) {
      infile>> vertexs[i];
    }
    for (auto i=0; i<b; i++) {
      infile>> positions[i];
    }
    for (auto i=0; i<c; i++) {
      infile>> dstAndPos[i];
    }
  }

  void gen_hybrid(delta_log<vertex> &log, graph<vertex> &graph, bool origin = true) {
    // cout << "start gen hybrid" << endl;
    // print_mem("BellmanFord");
    bool need_hybrid = graph.is_hybrid();
    if (need_hybrid && origin) {
      HVertex::carray.start_append();
    }
    uintT pos = 0;
    uintT next;
    while (pos < log.size()) {
      next = log.get_next(pos);
      uintE start = log.deltaLog[pos].first;
      vertex * vtx = graph.getvertex(start);
      if (vtx->is_high_degree_vertex()) {
        vertexs.push_back(start);
        positions.push_back(dstAndPos.size());
        bool waitformid = true;
        for (uintT i=pos; i<next; i++) {
          uintE data = log.deltaLog[i].second.first;
          intE weight = log.deltaLog[i].second.second;
          if (weight >= 0) {
            dstAndPos.push_back(data);
            dstAndPos.push_back(weight);
          } else {
            if (waitformid) {
              positions.push_back(dstAndPos.size());
              waitformid = false;
            }
            dstAndPos.push_back(data);
          }
        }
        if (waitformid) {
          positions.push_back(dstAndPos.size());
        }
      } else {
        // cout << "in low" << endl;
        vector<uintE> new_ver(vtx->getLastVersion());
        for (auto i=pos; i<next; i++) {
          uintE data = log.deltaLog[i].second.first;
          intE weight = log.deltaLog[i].second.second;
          if (weight >= 0) {
            vector_index_delete(new_ver, weight);
            // new_ver.index_delete(weight);
          } else {
            new_ver.push_back(data);
          }
        }
        vtx->append(new_ver, start, graph.get_version()+1);
      }
      pos = next;
    }
    if (need_hybrid && origin) {
      HVertex::carray.end_append();
    }
    // cout << "after gen delta : add" << get_add_size() << " and del " << get_del_number() << endl;
    // cout << "finish gen hybrid" << endl;
  }

  void write_edge_entry(const char * filename) {
    ofstream outfile;
    outfile.open(filename);
    if (!outfile.is_open()) {
      cout << "fail to open file " << filename << endl;
      abort();
    }  
    outfile << "DELTA_FILE" << endl;
    
    outfile << vs << " " << ve << endl;
    cout << get_add_size() << " " << get_del_number() << endl;
    outfile << get_add_size() << endl << get_del_number() << endl;

    for (auto i=0; i< vertexs.size(); i++) {
      auto from = vertexs[i];
      uintT mid, end;
      mid = positions[2*i+1];
      if (i == vertexs.size() - 1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      for (auto j=mid; j<end; j++) {
        outfile << from << " " << dstAndPos[j] << endl;
      }
    }
    for (auto i=0; i< vertexs.size(); i++) {
      auto from = vertexs[i];
      uintT start, mid;
      start = positions[2*i];
      mid = positions[2*i+1];
      for (auto j=start; j<mid; j+=2) {
        outfile << from << " " << dstAndPos[j] << endl;
      }
    }
    outfile.close();
  }

  void write_deltalog(const char * filename) {
    ofstream outfile;
    outfile.open(filename);
    if (!outfile.is_open()) {
      cout << "fail to open file " << filename << endl;
      abort();
    }
    outfile << "DELTA_LOG_FILE" << endl;
    outfile << vs << " " << ve << endl;
    outfile << get_add_size() + get_del_number() << endl;
    for (auto i=0; i< vertexs.size(); i++) {
      auto from = vertexs[i];
      uintT start, mid, end;
      start = positions[2*i];
      mid = positions[2*i+1];
      if (i == vertexs.size() - 1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      for (auto j=start; j<mid; j+=2) {
        outfile << from << " " << dstAndPos[j] << " " << dstAndPos[j+1] << endl;
      }
      for (auto j=mid; j<end; j++) {
        outfile << from << " " << dstAndPos[j] << " " << -1 << endl;
      }
    }
    outfile.close();
  }

  void write_pensieve_delta(const char * filename) {
    ofstream outfile;
    outfile.open(filename);
    if (!outfile.is_open()) {
      cout << "fail to open file " << filename << endl;
      abort(); 
    }
    cout << "begin to save file " << filename << endl;
    outfile << "PENSIEVE_DELTA_FILE" << endl;
    outfile << vs << " " << ve << endl;
    outfile << vertexs.size() << endl;
    outfile << positions.size() << endl;
    outfile << dstAndPos.size() << endl;
    for (auto i: vertexs) {
      outfile << i << endl;
    }
    for (auto i: positions) {
      outfile << i << endl;
    }
    for (auto i: dstAndPos) {
      outfile << i << endl;
    }
    outfile.close();
  }

  uintT get_path_end(vector<delta<vertex> *> &path, uintT ver_begin) {
    if (path.empty()) {
      cout << "error in get path end" << endl;
    }
    // cout << ver_begin << " to ";
    for (auto i: path) {
      if (i->vs == ver_begin) {
        ver_begin = i->ve;
      } else if (i->ve == ver_begin) {
        ver_begin = i->vs;
      }
      else {
        cout << "invalid path from " << endl;
        print_path(path);
      }
    }
    // cout << ver_begin << endl;
    return ver_begin;
  }

  // return all add pair of delta
  void pack_add_items(vector<pair<uintE, uintE>> * ret) {
    for (auto i=0; i<vertexs.size(); i++) {
      uintT start, end;
      start = positions[2*i+1];
      if (i == vertexs.size()-1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      if (end > dstAndPos.size()) {
        cout << i << " " << vs << " " << end << " " << start << " " << dstAndPos.size() << endl;
        abort();
      }
      for (auto j=start; j<end; j++) {
        ret->push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
  }

  // return all del pair of delta
  void pack_del_items(vector<pair<uintE, uintE>> * ret) {
    for (auto i=0; i<vertexs.size(); i++) {
      for (auto j=positions[i*2]; j<positions[i*2+1]; j+=2) {
        ret->push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
  }

  uintE get_delta_size() {
    return dstAndPos.size();
  }

  uintT get_del_number() {
    // if (del_size >= 0) {
    //   return (uintT) del_size;
    // }
    uintT ret=0;
    uintT num = positions.size() / 2;

    for (auto i=0; i<num; i++) {
      ret += (positions[2*i+1] - positions[2*i]) / 2;
    }
    del_size = ret;
    return ret;
  }

  uintT get_del_size() {
    return 2* get_del_number();
  }

  uintT get_add_size() {
    return dstAndPos.size() - get_del_size();
  }

  void print_vertex(uintT vertex_id, bool add=false) {
    for (auto i=0; i<vertexs.size(); i++) {
      if (vertexs[i] == vertex_id) {
        uintT s = positions[2*i];
        uintT m = positions[2*i+1];
        uintT e = i==vertexs.size()-1 ? dstAndPos.size() : positions[2*i+2];
        for (auto j=s; j<m; j+=2) {
          cout << dstAndPos[j] << " @ " << dstAndPos[j+1] << endl;
        }
        if (add) {
          for (auto j=m; j<e; j++) {
            cout << dstAndPos[j] << endl;
          }
        }
        break;
      }
    }
  }

  int find_edge(uintE s, uintE e) {
    for (auto i=0; i<vertexs.size(); i++) {
      if (vertexs[i] != s) {
        continue;
      }
      uintT begin, mid, end;
      begin = positions[2*i];
      mid = positions[2*i+1];
      if (i == vertexs.size() - 1) {
        end = dstAndPos.size();
      } else {
        end = dstAndPos[2*i+2];
      }
      for (auto j=begin; j<mid; j++) {
        if (dstAndPos[j] == e) {
          return dstAndPos[j+1];
        }
      }
      for (auto j=mid; j<end; j++) {
        if (dstAndPos[j] == e) {
          return -1;
        }
      }
    }
    return -2;
  }

  void shrink() {
    vertexs.shrink_to_fit();
    positions.shrink_to_fit();
    dstAndPos.shrink_to_fit();
  }

  void clear() {
    vector<uintE>().swap(vertexs);
    vector<uintE>().swap(positions);
    vector<uintE>().swap(dstAndPos);
  }
};

template<class vertex>
void print_path(vector<delta<vertex> *> path) {
  if (path.empty()) {
    cout << "empty path" << endl;
    abort();
    return;
  }
  for (auto i: path) {
    cout << i-> vs << " -> " << i->ve << endl;
  }
  cout << endl;
}

template <class vertex> 
vector<uintT> get_trace(vector<delta<vertex>> &deltas, uintT target) {
  vector<uintT> ret;
  while (target) {
    for (auto i=0; i<deltas.size(); i++) {
      if (deltas[i].ve == target) {
        ret.push_back(i);
        target = deltas[i].vs;
      }
    }
  }
  return ret;
}

template <class vertex> 
vector <uintT> get_revert_path(vector<delta<vertex>> &deltas, uintT from, uintT to) {
  if (from == to) {
    vector<uintT> ret;
    return ret;
  }

  vector<uintT> from_trace = get_trace(deltas, from);
  vector<uintT> to_trace = get_trace(deltas, to);

  while (!from_trace.empty() && !to_trace.empty() 
          && (from_trace[from_trace.size()-1] == (to_trace[to_trace.size()-1]))) {
    from_trace.pop_back();
    to_trace.pop_back();
  }

  reverse(to_trace.begin(), to_trace.end());

  from_trace.insert(from_trace.end(), to_trace.begin(), to_trace.end());
  return from_trace;
}

template <class vertex>
int apply(graph<vertex> & ga, delta<vertex> & da) {
  // version check
  if (ga.get_version() != da.vs) {
    std::cout << "version mismatch, delta apply failed, nothing changed" << endl;
    cout << "version = " << ga.get_version() << " and delta " << da.vs << " " << da.ve << endl;
    return -1;
  }
  // cout << "apply from " << da.vs << " to " << da.ve << endl;
  // int count_copy;
  // int count_add, count_del;
  // double time_copy, time_delta;
  // ChronoTimer cter;
  if (ga.is_hybrid()) {
    auto start = HVertex::carray.get_start_offset_of_version(ga.get_version());
    auto end = HVertex::carray.get_start_offset_of_version(da.ve);
    // count_copy = end - start;
    for (auto i=start; i<end; i++) {
      ga.getvertex(HVertex::carray.vertexs[i])->forward();
    }
  }
  // time_copy = cter.elapsed();
  
  int count = (int) da.vertexs.size();
  // cter.start();
  
  // print_mem("BellmanFord");
  for(int i=0; i<count; i++) {
    // get target vertex
    auto vtmp = ga.getvertex(da.vertexs[i]);

    // do the delete first
    for(int j=da.positions[2*i]; j<da.positions[2*i+1]; j+= 2) {
      vtmp->index_delete(da.dstAndPos[j+1]);
    }

    vector<uintE>::iterator it = da.dstAndPos.begin();
    std::size_t b, e;
      // then do some add
    if (i == count - 1) {
      b = da.positions[2*i+1];
      vtmp->push_back(it + b, da.dstAndPos.end());
    } else {
      b = da.positions[2*i+1];
      e = da.positions[2*i+2];
      vtmp->push_back(it+b, it+e);
    }
  }
  // time_delta = cter.elapsed();
  // cout << "apply finish" << endl;
  // cout << "copy_size " << count_copy << " with time " << time_copy << endl;
  // cout << "copy : " << time_copy/count_copy << endl;
  // cout << "delta size " << da.get_add_size() << " with time " << time_delta << endl;
  // cout << "delta : " << time_delta / (da.get_add_size()) << endl;
  
  ga.set_version(da.ve);
  // graph.update_m();
  int add = ((int)(da.get_add_size()) - (int)(da.get_del_number()));
  ga.add_m(add);
  return 0;
}

template <class vertex>
int revert(graph<vertex> &ga, delta<vertex> &da) {
  // version check
  if (ga.get_version() != da.ve) {
    std::cout << "version mismatch, delta revert failed, nothing changed" << endl;
    cout << "version = " << ga.get_version() << " and delta " << da.vs << " " << da.ve << endl;
    return -1;
  }

  // int count_copy;
  // int count_add, count_del;
  // double time_copy, time_delta;

  // cout << "revert from " << da.vs << " to " << da.ve << endl;
  // cout << "dirty data before " << ga.accessAllEdges() << endl;
  // ChronoTimer cter;
  if (ga.is_hybrid()) {
      // high_degree_vertex dont reponse this call and wait for following delta information
    int start = (int) HVertex::carray.get_start_offset_of_version(da.vs);
    auto end = HVertex::carray.get_start_offset_of_version(ga.get_version());
    // cout << start << " " << end << " " << HVertex::carray.vertexs.size() << endl;
    // count_copy = end - start;
    for (int i=end-1; i>=start; i--) {
      ga.getvertex(HVertex::carray.vertexs[i])->backward();
    }
  }
  // time_copy = cter.elapsed();
  
  int count = (int) da.vertexs.size();
  // cter.start();
  // cout << "after ld " << cter.elapsed() << endl;
  {parallel_for(int i=0; i<count; i++) {
  // {for(int i=0; i<count; i++) {
    vertex* vtmp = ga.getvertex(da.vertexs[i]);
    if (i == count - 1) {
      for(int j=da.positions[2*i+1]; j<da.dstAndPos.size(); j+= 1) {
        vtmp->pop_back();
      }
    } else {
      for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
        vtmp->pop_back();
      }
    }
    for(int j=da.positions[2*i+1]-2; j>=(int) da.positions[2*i]; j-= 2) {
      vtmp->index_addtion(da.dstAndPos[j], da.dstAndPos[j+1]);
    }
  }}
  // time_delta = cter.elapsed();
  // cout << "revert finish" << endl;
  // cout << "copy_size " << count_copy << " with time " << time_copy << endl;
  // cout << "copy : " << time_copy/count_copy << endl;
  // cout << "delta size " << da.get_add_size()  << " with time " << time_delta << endl;
  // cout << "delta : " << time_delta / (da.get_add_size()) << endl;
  // cout << "after index add " << cter.elapsed() << endl;
  ga.set_version(da.vs);
  int add = ((int)(da.get_add_size()) - (int)(da.get_del_number()));
  ga.add_m(-1*add);
  return 0;
}

template <class vertex>
int jump(graph<vertex>& ga, vector<delta<vertex> *>&path) {
  // if (!path.empty())
    // print_path<vertex>(path);
  // cout << ga.get_version() << endl;
  ChronoTimer cter;
  for (auto i=0; i<path.size(); i++) {
    auto d = path[i];
    // cout << d->vs << " to " << d->ve << " with graph " << ga.get_version() << endl;
    if (d->vs == ga.get_version()) {
      apply(ga, *d);
    } else if (d->ve == ga.get_version()) {
      revert(ga, *d);
    } else {
      cout << "error in jump with path" << endl;
    }
    // cout << cter.elapsed() << endl;
  }
  return 0;
}

#endif
