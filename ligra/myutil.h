#ifndef MY_UTIL_H
#define MY_UTIL_H

#include <chrono>
#include <ctime>
#include <random>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#define DEBUG 1

#ifdef DEBUG
#define V_V_V_V_V_V_logstart         if(1){
#define V_V_V_V_V_V_logif(cond)      if(cond){
#else
#define V_V_V_V_V_V_logstart         {
#define V_V_V_V_V_V_logif(cond)      {
#endif
#define V_V_V_V_V_V_logend           }
#define V_V_V_V_V_V_logprepare       ;
#define V_V_V_V_V_V_logpreparefinish ;


/*     chrono defines       */
class ChronoTimer {
public:
  ChronoTimer() : beg_(clock_::now()) {}
  void start() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<nano_>(clock_::now() - beg_).count();
  }
  double elapsed_milli() const {
    return std::chrono::duration_cast<milli_>(clock_::now() - beg_).count();
  }
  double moment() const {
    auto tp = clock_::now();
    return std::chrono::duration_cast<nano_>(tp.time_since_epoch()).count();
  }
  double moment_milli() const {
    auto tp = clock_::now();
    return std::chrono::duration_cast<milli_>(tp.time_since_epoch()).count();
  }

  double elapsed_and_refresh() {
    auto ret = elapsed();
    start();
    return ret;
  }

  double elapsed_and_refresh_milli() {
    auto ret = elapsed_milli();
    start();
    return ret;
  }

private:
  using clock_ = std::chrono::high_resolution_clock;
  using nano_ = std::chrono::nanoseconds;
  using milli_ = std::chrono::milliseconds;
  std::chrono::time_point<clock_> beg_;
};

uintT get_new_from(uintT current_max, uintT bound = 5) {
  if (bound == 0) {
    return current_max;
  }
  uintT current_min;
  if (current_max >= bound) {
    current_min = current_max - bound;
  } else {
    current_min = 0;
  }
  // std::cout << current_max << " " << current_min << std::endl;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(current_min, current_max);
  auto x = dis(gen);
  // std::cout << x << std:: endl << std::endl;
  return x;
}



/* func for index delete and index addition.
 *
 * 
 * 
*/

// try not use namespace in this file.
// using namespace std;

template <class T>
void vector_index_delete(std::vector<T> & vec, uintT pos) {
  if (vec.size() > pos) {
    vec[pos] = vec.back();
  }
  if (vec.size() > 0) {
    vec.pop_back();
  }
}

template <typename T>
void vector_index_addition(std::vector<T>& vec, uintT pos, const T data) {
  auto s = vec.size();
  if (pos < s) {
    vec.push_back(vec[pos]);
    vec[pos] = data;
  } else {
    vec.push_back(data);
  }
}


typedef std::pair<uintE, pair<uintE, intE>> intTriple;


// reload of () for sort, all for lessthan.
struct logLT {
  // sort by : 1, 3, 2
  bool operator () (intTriple a, intTriple b) {
    if (a.first != b.first) {
      return a.first < b.first;
    } else if (a.second.second != b.second.second) {
      return a.second.second > b.second.second;
    } else {
      return a.second.first < b.second.first;
    }
  }
  // sort by : 1, 2
  bool operator () (std::pair<uintT, uintT> a, std::pair<uintT, uintT> b) {
    if (a.first != b.first) {
      return a.first < b.first;
    }
    return a.second < b.second;
  }
  // sort by : 1, 2, 3+4, 3
  bool operator () (std::pair<std::pair<uintT, uintT>, std::pair<uintT, uintT>> a, 
                    std::pair<std::pair<uintT, uintT>, std::pair<uintT, uintT>> b) {
    if (a.first.first != b.first.first) {
      return a.first.first < b.first.first;
    } else if (a.first.second != b.first.second) {
      return a.first.second < b.first.second;
    } else if (a.second.first + a.second.second != b.second.first + b.second.second) {
      return a.second.first + a.second.first < b.second.first + b.second.second;
    } else {
      return a.second.first < b.second.first;
    }
  }
};


// get a random float and compare it with given level.
// usually use in check occur with given posibility.
bool randomFloatBiggerThan(double level) {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> ran(0, 1.0);
  return ran(gen) >= level;
}


// two func for compare.
bool equal(std::pair<uintE, uintE> a, std::pair<uintE, uintE> b) {
  return a.first == b.first && a.second == b.second;
}

bool lessthan(std::pair<uintE, uintE> a, std::pair<uintE, uintE> b) {
  return a.first < b.first || (a.first == b.first && a.second < b.second);
}

// function split
std::vector<std::string> split(const std::string a, const char delim = ' ') {
  std::vector<std::string> ret;
  std::size_t last = 0;
  std::size_t index = a.find_first_of(delim);
  while (index != std::string::npos) {
    if (index != last) {
      ret.push_back(a.substr(last, index));
    }
    last = index + 1;
    index = a.find_first_of(delim, last);
  }
  if (last < a.size()) {
    ret.push_back(a.substr(last));
  }
  return ret;
}

// split and get last
std::string split_last(const std::string a, const char delim = ' ') {
  auto rets = split(a, delim);
  if (!rets.empty()) {
    return rets.back();
  }
  return std::string();
}

const std::vector<std::string> short_path{"journal", "dewiki", "flickr", "enwiki", "orkut", "twitter"};
const std::vector<std::string> long_path{
  "/public/home/tangwei/graphdata/livejournal/",
  "/public/home/tangwei/graphdata/link-dynamic-dewiki/",
  "/public/home/tangwei/graphdata/flickr/",
  "/public/home/tangwei/graphdata/enwiki/",
  "/public/home/tangwei/graphdata/orkut/",
  "/public/home/tangwei/graphdata/twitter/"
};
const std::vector<std::string> delta_types{"smalldelta", "delta", "smalltree", "deletedelta"};

// transform csr path to delta path prefix.
std::string get_prefix(std::string s) {
  for (auto i=0; i<short_path.size(); i++) {
    if (s.find(short_path[i]) != s.npos) {
      return long_path[i];
    }
  }
  std::cout << "error in get prefix of csr dataset :" << endl << s << endl;
  abort();
}

std::string get_basedir(const std::string s, double delta_rate = 2, double add_rate = 2, int delta_type = 0) {
  auto prefix = get_prefix(s);
  if (s.find("flickr") == s.npos) {
    return prefix + "smalldelta/";
  }
  if (delta_rate != 2 || add_rate != 2) {
    if (delta_rate == 2) {
      delta_rate = 0.001;
    }
    if (add_rate == 2) {
      add_rate = 0.9;
    }
    ostringstream str1, str2;
    str1 << delta_rate;
    string delta_rate_string = str1.str();
    str2 << add_rate;
    string add_rate_string = str2.str();
    return prefix + "d" + delta_rate_string + "-a" + add_rate_string + '/';
  }
  std::cout << "delta type " << delta_types[delta_type] << endl;
  return prefix + delta_types[delta_type] + '/';
}

#endif