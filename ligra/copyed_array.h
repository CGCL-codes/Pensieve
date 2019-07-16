#ifndef _COPYED_ARRAY_H
#define _COPYED_ARRAY_H

#include <vector>

#include "utils.h"

using std::vector;
using std::cout;
using std::endl;

struct copyed_array
{
  vector<uintT> versions{0};
  vector<int> froms;
  vector<uintT> offsets{0};
  vector<uintE> vertexs;
  vector<uintT> nexts;
  vector<uintE> edges;
  // when this struct push back data, the open_append should be true, by start_append.
  // then append finish, set to false and update version info.
  bool open_append{false};

  copyed_array() {}

  // return new pos for first append of a vertex to modify pos to get its next.
  int append(vector<uintE> & data, uintE vtx, int version, int prevpos) {
    if (!open_append) {
      cout << "error to append to a copyed_array with open_append false." << endl;
      cout << "append issue rejected." << endl;
      return -1;
    }
    // the ret pos is just the last pos of offsets, not the after last.
    int ret = offsets.size()-1;
    if (prevpos == -1) {
      froms.push_back(-1 * ret - 1);
    } else {
      froms.push_back(prevpos);
      nexts[prevpos] = ret;
    }
    vertexs.push_back(vtx);
    edges.insert(edges.end(), data.begin(), data.end());
    offsets.push_back(edges.size());
    nexts.push_back(-1);
    return ret;
  }

  uintT get_start_offset_of_version(uintT ver) {
    if (ver >= versions.size()) {
      cout << "error in get start offset of version " << ver 
           << " with " << versions.size() << " versions" << endl;
    return 0;
    }
    return versions[ver];
  }

  uintT get_end_offset_of_version(uintT ver) {
    if (ver >= versions.size()) {
      cout << "error in get end offset of version " << ver 
           << " with " << versions.size() << " versions" << endl;
      return 0;
    }
    return versions[ver+1];
  }

  void start_append() {
    if (!open_append) {
      open_append = true;
    } else {
      cout << "error to reopen a copyed_array to append" << endl;
    }
  }

  void end_append() {
    if (open_append) {
      versions.push_back(offsets.size()-1);
      open_append = false;
    } else {
      cout << "error to close a non-open copyed_array" << endl;
    }
  }
};


#endif