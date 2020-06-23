#include "delta.h"
#include "IO.h"

using namespace std;

int main(void) {
  graph<symmetricVertex> G = readGraph<symmetricVertex>(
      "", 
      0, 1, 0, 0);
  int count = 10;
  for(auto i = 0; i < count; i++) {
    delta_log<symmetricVertex> dlg = delta_log(G);
    delta<symmetricVertex> dlt = delta(dlg, G);
    writeLog(dlg, dir, count);
    writeDelta(dlt, dir, count);
    apply(G, dlt);
  }
  return 0;
}




/* #region *
graph g = load(filename) ///
log = gen_log(g, parameter)
delta = gen_delta(g, log)///
apply(graph, delta)///
bigdelta.append(delta)
delta.writedown()

* #endregion */