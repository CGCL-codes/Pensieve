#ifndef _EXP_DIST_H
#define _EXP_DIST_H

#include <random>
#include <ctime>
using namespace std;

const int DIST_SIZE=20;

struct expDist {
  static default_random_engine e;
  
  expDist() {
    e.seed(time(0) * time(0));
  }
  
//   float getRaw() {
//     return expd(e);
//   }
  int getRand(int maxrand) {
    int ret;
    exponential_distribution<float> expd(0.5);
    do {
      ret = (int) lround(10 * expd(e));
    }
    while (ret >= maxrand);
    return maxrand - ret;
  }
};

default_random_engine expDist::e;

struct distOrder {
  static int order[DIST_SIZE];

  static int get_size() {
    return DIST_SIZE;
  }

  static int get(int index, int want_max=60) {
    if (index >= DIST_SIZE || index < 0) {
      cout << "index out of range " << index << endl;
      return 0;
    }
    if (want_max == 60) {
      return order[index];
    }
    auto tmp = order[index];
    return int(order[index]*want_max / 60);
  }
};

int distOrder::order[DIST_SIZE]{60, 59, 54, 43, 55, 58, 55, 52, 58, 49, 
                                56, 56, 43, 49, 52, 58, 38, 45, 59, 52};

#endif
