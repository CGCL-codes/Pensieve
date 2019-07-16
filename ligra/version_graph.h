#ifndef _VERSION_GRAPH_H
#define _VERSION_GRAPH_H

#include <vector>
#include "graph.h"
#include "delta.h"

// for -g (version graph)
template <class vertex>
struct versionGraph {
  int divide = 5;
  vector<delta<vertex>*> tree;
  vector<delta<vertex>*> chain;
  vector<pair<uintT, uintT>> edges;

  versionGraph() {}

  versionGraph(float di) {
    if (di > 1) {
      divide = (int) di;
    }
  }

  vector<delta<vertex> *> get_trace(uintT target) {
    vector<delta<vertex> *> ret;
    bool edit = true;
    while (target && edit) {
      edit = false;
      for (auto i=0; i<edges.size(); i++) {
        if (edges[i].second == target) {
          ret.push_back(get_edge(edges[i].first, edges[i].second));
          target = edges[i].first;
          edit = true;
        }
      }
    }
    return ret;
  }

  vector <delta<vertex> *> get_revert_path(uintT from, uintT to) {
    if (from == to) {
      vector<delta<vertex> *> ret;
      return ret;
    }

    vector<delta<vertex> *> from_trace = get_trace(from);
    vector<delta<vertex> *> to_trace = get_trace(to);

    while (!from_trace.empty() && !to_trace.empty() 
          && (from_trace[from_trace.size()-1]->vs == to_trace[to_trace.size()-1]->vs)
          && (from_trace[from_trace.size()-1]->ve == to_trace[to_trace.size()-1]->ve)) {
      from_trace.pop_back();
      to_trace.pop_back();
    }

    reverse(to_trace.begin(), to_trace.end());

    from_trace.insert(from_trace.end(), to_trace.begin(), to_trace.end());
    vector<delta<vertex> *> tmp;
    for (auto i: from_trace) {
      tmp.push_back(i);
    }
    return tmp;
  }

  // update tree, often just pushback, sometimes add a totally new delta and update chain
  void update_tree(graph<vertex> &graph, delta<vertex> * moved) {
    // try to get a vector contains distance between new version and others. 

    ChronoTimer cter;
    auto current_upper = get_new_upper(moved);
    
    if (current_upper == moved->vs) {
      tree.push_back(moved);
      chain.erase(chain.begin());
    } else {
      auto path = get_revert_path(graph.get_version(), current_upper);
      jump(graph, path);
      
      path = get_revert_path(current_upper, moved->ve);
      
      delta<vertex> *new_delta = new delta<vertex> (graph, path, current_upper, moved->ve);

      edges.push_back(make_pair(current_upper, moved->ve));
      remove_edge(moved->vs, moved->ve);
      cout << "new edge " << current_upper << " => " << moved->ve
           << " instead of " << moved->vs << " -> " << moved->ve << endl;
      delete moved;
      
      tree.push_back(new_delta);
      
      
      apply(graph, *new_delta);

      update_chain(graph);
    }
  }

  void remove_edge(uintT s, uintT e) {
    for (auto i=edges.begin(); i != edges.end(); i++) {
      if ((i->first == s && i->second == e) || (i->first == e && i->second == s)) {
        edges.erase(i);
        return;
      }
    }
    cout << "fail to erase edge from " << s << " to " << e << endl;
  }

  // update only if new delta can not be simply pushback 
  // loop all the chain and reformat delta 
  void update_chain(graph<vertex> &graph) {
    if (graph.get_version() != chain[1]->vs) {
      cout << "error in update chain, graph version is " << graph.get_version() << " and chain from " << chain[0]->vs << endl; 
      abort();
    }
    vector<delta<vertex> *> new_chain;

    // chain[0] is previous deleted in update_tree
    for (auto i=1; i<chain.size(); i++) {
      vector<delta<vertex> *> tmp;
      tmp.push_back(chain[i]);
      auto new_delta = new delta<vertex>(graph, tmp, graph.get_version(), graph.get_version() + 1);

      new_chain.push_back(new_delta);
      delete chain[i];
      apply(graph, *new_delta);
    }
    chain = new_chain;
    // print_edges();
  }

  // get max version end of all version graph(include tree and chain)
  // version of delta of chain is always front of tree.
  // 获取整个版本图的最大版本。
  uintT get_max_version() {
    vector<uintT> ver_ends;
    for (auto i : chain) {
      ver_ends.push_back(i->ve);
    }
    return *max_element(ver_ends.begin(), ver_ends.end());
  }

  // get neighbors of a version in tree.
  // todo: diff tree and chain
  vector<uintT> get_ngh(uintT ver) {
    vector<uintT> ret;
    for (auto edge: edges) {
      if (edge.first == ver) {
        ret.push_back(edge.second);
      } else if (edge.second == ver) {
        ret.push_back(edge.first);
      }
    }
    return ret;
  }

  // get not accurate delta size
  int estimate_delta(vector<delta<vertex> *> path) {
    int ret = 0;
    for (auto i: path) {
      ret += i->get_add_size();
      ret -= i->get_del_size();
    }
    return ret > 0 ? ret : ret * -2;
  }

  // get accurate delta size
  int accurate_delta(vector<delta<vertex> *> path) {
    vector<pair<uintE, uintE>> ret_add, ret_del;
    for (auto i : path) {
      i->pack_add_items(&ret_add);
      i->pack_del_items(&ret_del);
    }
    sort(ret_add.begin(), ret_add.end(), logLT());
    sort(ret_del.begin(), ret_del.end(), logLT());
    // cout << ret_add.size() << " " << ret_del.size() << endl;
    uintT pa = 0, pd = 0;
    uintT count_add=0, count_del=0;
    while (pa < ret_add.size() && pd < ret_del.size()) {
      if (equal(ret_add[pa], ret_del[pd])) {
        pa++; pd++;
      } else if (lessthan(ret_add[pa], ret_del[pd])) {
        count_add++;
        pa++;
      } else {
        pd ++;
        count_del++;
      }
    }
    if (pa < ret_add.size()) {
      count_add += ret_add.size() - pa;
    } else if (pd < ret_del.size()) {
      count_del += ret_del.size() - pd;
    } else if (pa > ret_add.size() || pd > ret_del.size()){
      cout << "something error in accurate_del" << endl;
    }
    
    return count_add + count_del * 2;
  }

  // get where the new edge of graph insert. 
  uintT get_new_upper(delta<vertex> * moved) {
    uintT tree_size = tree.size();

    int current_min = moved->get_delta_size();
    
    uintT current_upper = tree_size;

    // moved.vs == tree_size  no need to loop from tree_size
    for (int i=(int)tree_size-1; i>=0; i--) {
      
      vector<delta<vertex> *> path = get_revert_path(moved->vs, i);
      if (path.empty()) {
        continue;
      }
      path.push_back(moved);
      auto tmp = estimate_delta(path);
      // print_path<vertex>(path);
      // cout << "estimate " << moved->vs << " -> " << i << " " << tmp << " " << current_min << endl;
     
      if (tmp < current_min) {
        int new_min = accurate_delta(path);
        // cout << "accurate " << moved->vs << " -> " << i << " " << new_min << " " << current_min << endl;
        if (new_min < current_min) {
          current_min = new_min;
          current_upper = i;
        }
      }
    }
    return current_upper;
  }

  // get delta by version start and version end; return null if not found
  delta<vertex> * get_edge(uintT b, uintT e) {
    if (b == e) {
      cout << "get error edge from " << b << " to " << e << endl;
      print_edges();
      abort();
      return NULL;
    }

    for (auto i = 0; i < tree.size(); i++) {
      auto t = tree[i];
      if ((t->vs == b && t->ve == e)) {
        return t;
      }
    }
    for (auto i = 0; i < chain.size(); i++) {
      auto c = chain[i];
      if ((c->vs == b && c->ve == e)) {
        return c;
      }
    }
    cout << "not found edge from " << b << " to " << e << endl;
    print_edges();
    return NULL;
  }

  void append(graph<vertex> & graph, delta<vertex>* d) {
    chain.push_back(d);
    edges.push_back(make_pair(d->vs, d->ve));
    repart(graph);
    // print_edges();
  }

  void repart(graph<vertex> & graph) {
    if (chain.size() > divide) {
      delta<vertex> * moved = chain[0];
      update_tree(graph, moved);
    } else {
      cout << "chain not reach limit  " << divide << endl;
    }
  }

  void print_edges() {
    if (edges.empty()) {
      cout << "empty version graph " << endl;
      return;
    }
    cout << "edges:" << endl;
    for (auto i:edges) {
      cout << i.first << " -> " << i.second << endl;
    }
    cout << "tree:" << endl;
    for (auto i:tree) {
      cout << i->vs << " -> " << i->ve << endl;
      cout << i->get_add_size() << " vs " << i->get_del_number() << endl;
    }
    cout << "chain:" << endl;
    for (auto i:chain) {
      cout << i->vs << " -> " << i->ve << endl;
      cout << i->get_add_size() << " vs " << i->get_del_number() << endl;
    }
  }

  void shrink() {
    for (auto i: tree) {
      i->shrink();
    }
    for (auto i: chain) {
      i->shrink();
    }
  }
};

template <class vertex>
int jump(graph<vertex>& ga, versionGraph<vertex> & vg, int target) {
  // first update for initial version.
  // vg.update_vertex_size(ga);
  // cout << ga.get_version() << " " << target << endl;
  if (ga.get_version() == target) {
    cout << "jump to same version, nothing happened" << endl;
    return 0;
  } 
  auto path = vg.get_revert_path(ga.get_version(), target);
  // cout << "before print path " << endl;
  // print_path<vertex>(path);
  // cout << "after print path " << endl;
  if (path.empty()) {
    cout << "fail to jump from " << ga.get_version() << " to " << target <<endl;
  }
  jump(ga, path);
  // vg.update_vertex_size(ga);
  return 0;
}

#endif