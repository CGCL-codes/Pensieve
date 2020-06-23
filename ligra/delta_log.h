#ifndef _DELTA_LOG_H
#define _DELTA_LOG_H

#include <vector>
#include <unordered_set>
#include "vertex.h"
#include "myutil.h"
#include "threadpool.h"

using namespace std;

uintT get_next(vector<uintT> &arr, uintT target, uintT near) {
  uintT m = arr.size();
  if (near == 0) {
    near = 1;
  }
  while (arr[near] == target) {
    if (near == m-1) {
      return near;
    }
    near ++;
  }
  return near-1;
}

uintT bin_search(vector<uintT> &arr, uintT target) {
  uintT start = 0;
  uintT end = arr.size() - 1;
  while(end - start > 1) {
    if (arr[end] == target) {
      return get_next(arr, target, end);
    }
    if (arr[start] == target) {
      return get_next(arr, target, start);
    }
    uintT middle = (end+start)/2;
    if (arr[middle] > target) {
      end = middle;
    } else {
      start = middle;
    }
  }
  return start;
}

mutex mtx;

template <class vertex>
struct delta_log{
  vector<intTriple> deltaLog;
  uintT ver;
  uintT ver_end;
  double add_rate = 0.1;
  double delta_rate = 0.001;
  

  delta_log(vector<intTriple> _deltalog):deltaLog(_deltalog){}

  delta_log() {}

  delta_log(graph<vertex> &graph, double _add_rate, double _delta_rate, bool diff) {
    add_rate = _add_rate;
    delta_rate = _delta_rate;
    do_gen_deltalog(graph);
  }

  delta_log(graph<vertex> &graph) {
    do_gen_deltalog(graph);
  }

  delta_log(graph<vertex> &graph, double _add_rate, double _delta_rate) {
    add_rate = _add_rate;
    delta_rate = _delta_rate;
    do_gen_deltalog(graph);
  }

  void do_gen_deltalog(graph<vertex> &graph) {
    ver = graph.get_version();
    ver_end = ver + 1;
    int a=0, d=0;
    for(uintT i = 0; i < graph.n; i++) {
      // if (!(i%(graph.n/10))) {
      //   cout << "tick" << endl;
      // }
      auto vtmp = graph.getvertex(i);
      auto vsize = vtmp->getInDegree();
      auto add_and_del = get_delta_easy(vtmp, graph.n);

      for (auto deledge : add_and_del) {
        deltaLog.push_back(make_pair(i, make_pair(deledge.first, deledge.second)));
        if (deledge.second >= 0) {
          d += 1;
        } else {
          a += 1;
        }
      }
    }
    sort(deltaLog.begin(), deltaLog.end(), logLT());
    // cout << "when do gen deltalog : add " << a << " and del " << d << endl;
  }

  vector<pair<uintE, intE>> get_delta_easy(vertex* v, uintE max_id) {
    vector<pair<uintE, intE>> ret;
    uintT s = v->getInDegree();
    if (s == 0) {
      return ret;
    }
    int remain = s*delta_rate + 0.5;
    int del_number = remain*(1-add_rate)+0.5;
    int add_number = remain - del_number;
    unordered_set<int> nghs;
    for (auto i=0; i<s; i++) {
      nghs.insert(v->getInNeighbor(i));
    }
    for (auto i=rand()%(max_id/2); add_number > 0&&i<max_id; i++) {
      if (nghs.find(i) == nghs.end()) {
        add_number -= 1;
        ret.push_back({i, -1});
      }
    }
    for (auto i=0; i<del_number; i++) {
      ret.push_back({v->getInNeighbor(i), i});
    }
    return ret;
  }

  vector<pair<uintE, intE>> get_delta(vertex * v, uintE max) {
    vector<pair<uintE, intE>> ret;
    uintT s = v->getInDegree();
    // if there is no edge in the vertex, there is a possibility to get a new one.
    if (s == 0) {
      if (randomFloatBiggerThan(1-100.0/max)) {
        ret.push_back(make_pair(rand()%max, -1));
      }
      return ret;
    }

    // delta entry remained.
    double remain = s * delta_rate;
    vector<uintT> add_new, del_new;
    unordered_set<int> nghs;
    if (s>100) {
      for (auto i=0; i<s; i++) {
        nghs.insert(v->getInNeighbor(i));
      }
    }

    while (remain > 0) {
      // if remain less than 1, there is a possibility to get a new delta.
      if (remain < 1 && randomFloatBiggerThan(remain)) break;
      // if try to get a delta entry 
      if (randomFloatBiggerThan(1-add_rate)) {
        uintT tmp = rand() % max;
        if (s>100 && nghs.find(tmp) != nghs.end()){
          continue;
        }
        // if selected end vertex is already in the vertex, continue to find another one.
        else if (s < 100 && v->find_val(tmp) != -1) {
          continue;
        } else {
          add_new.push_back(tmp);
        }
      } else {
        // del_new gathers pos of edges try to delete.
        del_new.push_back(rand()%s);
      }
      remain -= 1;
    }
    // sort and erase for dedup.
    sort(add_new.begin(), add_new.end());
    add_new.erase(unique(add_new.begin(), add_new.end()), add_new.end());
    sort(del_new.begin(), del_new.end());
    del_new.erase(unique(del_new.begin(), del_new.end()), del_new.end());
    for (auto i: add_new) {
      ret.push_back(make_pair(i, -1));
    }
    for (auto i: del_new) {
      ret.push_back(make_pair(v->getInNeighbor(i), i));
    }
    return ret;
  }

  // get delta_log from exist add and delete delta entry without pos.
  delta_log(graph<vertex> & graph, vector<pair<uintE, uintE>> &add, vector<pair<uintE, uintE>> &del) {
    for (auto i : add) {
      deltaLog.push_back(make_pair(i.first, make_pair(i.second, -1)));
    }
    for (auto i : del) {
      intE pos = graph.getvertex(i.first)->find_val(i.second);
      deltaLog.push_back(make_pair(i.first, make_pair(i.second, pos)));
    }
    delta_rate = 1.0*deltaLog.size() / graph.get_edge_number();
    add_rate = 1.0 * add.size() / deltaLog.size();
    sort(deltaLog.begin(), deltaLog.end(), logLT());
  }

  int size() {
    return deltaLog.size();
  }

  // get next index while the from vertex unchange.
  uintT get_next(uintT current) {
    for (uintT i=current+1; i<size(); i++) {
      if (deltaLog[i].first != deltaLog[current].first) {
        return i;
      }
    }
    return size();
  }

  
/*  file structure
 *  1: DELTA_LOG_FILE
 *  2: (base version and target version)
 *  3: (size of deltalog)
 *  4: (first int triple of deltalog)
 *  ......
 *  end of file
 * */

  int write_deltalog_to_file(const char * filename) {
    return write_deltalog_to_file(string(filename));
  }

  int write_deltalog_to_file(const string & filename) {
    ofstream outfile;
    outfile.open(filename);
    if (!outfile.is_open()) {
      cout << "fail to open file " << filename<< endl;
      // abort();
      return -1;
    }
    outfile << "DELTA_LOG_FILE" << endl;
    outfile << ver << " " << ver_end << endl;
    outfile << deltaLog.size() << endl;
    for (auto const dl : deltaLog) {
      outfile << dl.first << " " << dl.second.first << " " << dl.second.second << endl;
    }
    outfile.close();
    return 0;
  }
};

template <class vertex>
delta_log<vertex> load_deltalog_from_file(graph<vertex> & g, const char * filename) {
  ifstream infile;
  infile.open(filename);
  if (!infile.is_open()) {
    cout << "error in open file " << filename << endl;
    abort();
  }
  string code, version_info;

  getline(infile, code);
  delta_log<vertex> d;
  // infile >> version_info;
  getline(infile, version_info);
  // there may be only version_start or both version_start and version_end.
  auto split = version_info.find(" ");
  if (split != version_info.npos) {
    d.ver_end = atoi(version_info.substr(split).c_str());
    d.ver = atoi(version_info.substr(0, split).c_str());
  } else {
    d.ver = atoi(version_info.c_str());
    d.ver_end = d.ver + 1;
  }
  uintE start, end;
  int pos;
  if (!code.compare("DELTA_LOG_FILE")) {
    uintT length;
    infile >> length;
    for (auto i=0; i < length; i++) {
      infile >> start >> end >> pos;
      d.deltaLog.push_back(make_pair(start, make_pair(end, pos)));
    }
  } else if (!code.compare("DELTA_FILE")) {
    vector<pair<uintE, uintE>> add;
    vector<pair<uintE, uintE>> del;
    uintT add_length, del_length;
    infile >> add_length;
    infile >> del_length;

    for (size_t i=0; i < add_length; i++) {
      infile >> start >> end;
      d.deltaLog.push_back(make_pair(start, make_pair(end, -1)));
    }
    for (size_t i=0; i < del_length; i++) {
      infile >> start >> end;
      pos = g.getvertex(start)->find_val(end);
      if (pos == -1) {
        // cout << "not found edge " << end << " in vertex " << start << " 's neighbor in version" << d.ver << endl; 
        // abort();
        // continue;
      }
      d.deltaLog.push_back(make_pair(start, make_pair(end, pos)));
    }
    quickSort(d.deltaLog.data(), d.deltaLog.size(), logLT());
    d.add_rate = add_length / (add_length + del_length);
    d.delta_rate = (add_length + del_length) / g.m;
  } else {
    cout << "bad delta file" << filename << " with code" << code << endl;
    abort();
  }
  infile.close();
  return d;
}

#endif
