#ifndef _BIG_DELTA_H
#define _BIG_DELTA_H

#include <vector>
#include "graph.h"

// for no -g
template <class vertex>
struct bigDelta {
  vector<uintE> versions;
  vector<uintE> vertexs;
  vector<uintE> positions;
  vector<uintE> dstAndPos;

  vector<int> add_number;

  int append(delta<vertex> &da) {
    versions.push_back(vertexs.size());
    vertexs.reserve(vertexs.size()+da.vertexs.size());
    
    for (auto i = 0; i < da.vertexs.size(); i++) {
      vertexs.push_back(da.vertexs[i]);
    }
    
    positions.reserve(positions.size() + da.positions.size());
    int last_length = dstAndPos.size();
    
    for (auto i = 0; i < da.positions.size(); i++) {
      positions.push_back(da.positions[i]+last_length);
    }

    dstAndPos.reserve(last_length + da.dstAndPos.size());
    for (auto i = 0; i < da.dstAndPos.size(); i++) {
      dstAndPos.push_back(da.dstAndPos[i]);
    }

    int b = (int)da.get_add_size();
    int c = (int)da.get_del_number();
    int a =  b-c ;

    add_number.push_back(a);

    return 0;
  }

  int get_max_version() {
    return versions.size();
  }

  uintT size() {
    return dstAndPos.size();
  }

  uintT capacity() {
    return dstAndPos.capacity();
  }

  int get_add_number(uintT ver) {
    if (ver >= add_number.size()) {
      cout << "ver " << ver << " size " << add_number.size() << endl;
      abort();
    }
    if (add_number.size() != versions.size()) {
      cout << "add " << add_number.size() << " version " << versions.size() << endl;
      abort();
    }
    return add_number[ver];
  }

  int get_add_number(uintT start, uintT end) {
    if (start == end) {
      return 0;
    } else if (start > end) {
      return get_add_number(end, start) * -1;
    }
    
    int ret = 0;
    for (auto i=start; i<end; i++) {
      ret += get_add_number(i);
    }
    
    return ret;
  }
};


template <class vertex> 
int forward(graph<vertex> &ga, bigDelta<vertex> &bda, int step = 1) {
  if (step == 0) {
    // cout << "forward 0 step, that's to say, nothing happened" << endl;
    return 0;
  }else if (step < 0) {
    return backward(ga, bda, -1 * step);
  }

  int version_start = ga.get_version();
  int version_end = version_start + step;
  int version_max = bda.get_max_version();
  if (version_end > version_max + 1) {
    cout << "try to get far version than exist, nothing happened." << endl;
    cout << version_end << " > " << bda.get_max_version() << endl;
    return -1;
  }

  if (ga.is_hybrid()) {
    auto start = HVertex::carray.get_start_offset_of_version(version_start);
    auto end = HVertex::carray.get_start_offset_of_version(version_end);
    for (auto i=start; i<end; i++) {
      auto vtx = HVertex::carray.vertexs[i];
      ga.getvertex(vtx)->forward();
    }
  }

  uintT add_count = 0, del_count = 0;
  // cout << version_start << " to " << version_end << endl;
  for (auto ver = version_start; ver < version_end; ver++) {
    int count = (ver == version_max - 1 ? 
                  bda.vertexs.size() - bda.versions[ver] :
                  bda.versions[ver+1] - bda.versions[ver]);

    int prefix = bda.versions[ver];

    {parallel_for(int i=prefix; i<prefix + count; i++) {
      vertex* vtmp = ga.getvertex(bda.vertexs[i]);

      // do the delete first
      for(int j=bda.positions[2*i]; j<bda.positions[2*i+1]; j+= 2) {
        vtmp->index_delete(bda.dstAndPos[j+1]);
        add_count ++;
      }
      vector<uintE>::iterator it = bda.dstAndPos.begin();
      std::size_t b, e;
      // then do some add
      if (ver == version_max - 1 && i == prefix + count - 1) {
        b = bda.positions[2*i+1];
        vtmp->push_back(it + b, bda.dstAndPos.end());
        del_count += bda.dstAndPos.size() - b;
      } else {
        b = bda.positions[2*i+1];
        e = bda.positions[2*i+2];
        vtmp->push_back(it+b, it+e);
        del_count += e-b;
      }
    }}
  }
  ga.set_version(version_end);
  // ga.update_m();
  ga.add_m(bda.get_add_number(version_start, version_end));
  return 0;
}

template <class vertex> 
int backward(graph<vertex> &ga, bigDelta<vertex> &bda, int step) {
  if (step == 0) {
    cout << "backward 0 step, that's to say, nothing happened" << endl;
    return 0;
  }else if (step < 0) {
    return forward(ga, bda, -1 * step);
  }
  
  int version_start = ga.get_version();
  if (version_start > bda.get_max_version() + 1) {
    cout << "try to travel from some version not exists: " << bda.get_max_version() << endl;
    return -1;
  }
  int version_end = version_start - step;
  if (version_end < 0) {
    cout << "try to get negative version " << version_end << ", nothing happened." << endl;
    return -1;
  }
  int version_max = bda.get_max_version();
  // cout << version_start << " " << version_end << endl;
  // cout << ga.is_hybrid() << endl;
  ChronoTimer cter;
  if (ga.is_hybrid()) {
    auto start = HVertex::carray.get_start_offset_of_version(version_start);
    auto end = HVertex::carray.get_end_offset_of_version(version_end);
    for (int i=end-1; i>=start; i--) {
      ga.getvertex(HVertex::carray.vertexs[i])->backward();
    }
  }
  // cout << "after ld " << cter.elapsed() << endl;
  // cout << version_start << " " << version_end << endl;
  
  for (int ver = version_start-1; ver > version_end-1; ver--) {
    int count = (ver == version_max - 1) ? 
                  bda.vertexs.size() - bda.versions[ver] :
                  bda.versions[ver+1] - bda.versions[ver];

    int prefix = bda.versions[ver];

    {parallel_for(int i=prefix; i<count+prefix; i++) {
      vertex* vtmp = ga.getvertex(bda.vertexs[i]);
      if (ver == version_max - 1 && i == prefix+count-1) {
        for(int j=bda.positions[2*i+1]; j<bda.dstAndPos.size(); j+= 1) {
          vtmp->pop_back();
        }
      }else {
        for(int j=bda.positions[2*i+1]; j<bda.positions[2*i+2]; j+= 1) {
          vtmp->pop_back();
        }
      }
      for(int j=bda.positions[2*i+1]-2; j >= (int) bda.positions[2*i]; j-= 2) {
        vtmp->index_addtion(bda.dstAndPos[j], bda.dstAndPos[j+1]);
      }
    }}
  }
  
  // cout << "after add and delete" << cter.elapsed() << endl;
  ga.set_version(version_end);

  // ga.update_m();
  int diff_m = bda.get_add_number(version_start, version_end);

  ga.add_m(diff_m);
  // cout << "after add m" << cter.elapsed() << endl;
  return 0;
}

template <class vertex>
int jump(graph<vertex>& ga, bigDelta<vertex>& bda, int target) {
  if (target < 0 || target > bda.get_max_version()) {
    cout << "try to jump to version not exist, nothing happened" << endl;
    cout << target << "   " << bda.get_max_version() << endl;
    return -1;
  }
  // cout << "jump " << target << " " << ga.get_version() << endl;
  return forward(ga, bda, target - (int) ga.get_version());
}

template <class vertex>
class deltavector{
  vector<int> in_memory;
  vector<delta<vertex>> deltas;
  string savedir;

public:

  deltavector(string save = "/public/home/tangwei/graphdata/twitter/disktmp") {
    savedir = save;
  }
  void append(delta<vertex>& d) {
    deltas.push_back(d);
    // cout << "copy construction" << endl;
    in_memory.push_back(1);
  }

  void append(delta<vertex>&& d) {
    deltas.push_back(d);
    cout << "move construction" << endl;
    in_memory.push_back(1);
  }

  void save_disk(int count = 0) {
    if (deltas.size() < 10) {
      return;
    }
    if (count == 0) {
      count = deltas.size() * 0.8;
    }
    cout << "save from 0 to " << count << endl;;
    for (auto i=0; i<count; i++) {
      if (!in_memory[i]) {
        continue;
      }
      string save_file = get_disk_filename(i);
      deltas[i].write_pensieve_delta(save_file.c_str());
      deltas[i].clear();
      in_memory[i] = 0;
      cout << "saved " << save_file << endl;
    }
  }

  int save_one(bool trick = false) {
    auto remain = count_if(in_memory.begin(), in_memory.end(), [](int i){return i==1;});
    if (remain < 10 || remain <= in_memory.size() * 0.8) {
      return -1;
    }
    auto tmp = find(in_memory.begin(), in_memory.end(), 1);
    if (tmp == in_memory.end()) {
      return -1;
    }
    auto index = *tmp;
    if (!trick) {
      string save_file = get_disk_filename(index);
      deltas[index].write_pensieve_delta(save_file.c_str());
    }
    deltas[index].clear();
    in_memory[index] = 1;
    return 0;
  }

  void trick_save_disk(int count = 0) {
    if (deltas.size() < 10) {
      return;
    }
    if (count == 0) {
      count = deltas.size() * 0.8;
    }
    cout << "fake save from 0 to " << count << endl;
    for (auto i=0; i<count; i++) {
      if (!in_memory[i]) {
        continue;
      }
      deltas[i].clear();
      in_memory[i] = 0;
    }
  }

  bool delta_in_memory(uintT index) {
    return index < in_memory.size() && in_memory[index] == 1;
  }

  delta<vertex> & get_delta(uintT index) {
    return deltas[index];
  }

  string get_disk_filename(uintT index) {
    return savedir + "/" + to_string(index) + "-" + to_string(index+1);
  }
};

template <class vertex>
void jump_with_disk_check(deltavector<vertex> & dv, graph<vertex> & g, int target) {
  auto curr_ver = g.get_version();
  cout << "jump from " << curr_ver << " to " << target << endl;
  if (target == curr_ver) return;
  if (target>curr_ver) {
    for (auto i=curr_ver; i<target; i++) {
      if (dv.delta_in_memory(i)) {
        cout << "apply in memory " << i << endl;
        apply(g, dv.get_delta(i));
      } else {
        cout << "apply disk " << i << endl;
        delta<vertex> tmpd(dv.get_disk_filename(i).c_str());
        apply(g, tmpd);
      }
    }
  } else {
    for (int i=curr_ver-1; i>=target; i--) {
      if (dv.delta_in_memory(i)) {
        cout << "revert in memory " << i << endl;
        revert(g, dv.get_delta(i));
      } else {
        cout << "revert disk " << i << endl;
        delta<vertex> tmpd(dv.get_disk_filename(i).c_str());
        revert(g, tmpd);
      }
    }
  }
}

#endif
