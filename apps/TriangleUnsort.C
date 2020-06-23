#include "ligra.h"

template <class vertex>
long countCommon(vertex& A, vertex& B, uintE a, uintE b) {
    uintT nA = A.getInDegree();
    uintT nB = B.getInDegree();
    uintE * nghA = (uintE*)A.getInNeighbors();
    uintE * nghB = (uintE*)B.getInNeighbors();
    long ans = 0;
    
    if (nA < 100 && nB < 100) {
        for (auto i=0; i<nA; i++) {
            if (nghA[i]>=a) continue;
            for (auto j=0; j<nB; j++) {
                if (nghB[j] == nghA[i] && nghB[j]<b) {
                    ans ++;
                    break;
                }
            }
        }
    } else if (nA>=nB) {
        unordered_set<int> setA;
        for (auto i=0; i<nA; i++) {
            if (nghA[i]<a) {
                setA.insert(nghA[i]);
            }
        }
        for (auto i=0; i<nB; i++) {
            if (nghB[i] < b && setA.find(nghB[i])!=setA.end()) {
                ans += 1;
            }
        }
    } else {
        unordered_set<int> setB;
        for (auto i=0; i<nB; i++) {
            if (nghB[i]<b) {
                setB.insert(nghB[i]);
            }
        }
        for (auto i=0; i<nA; i++) {
            if (nghA[i] < a && setB.find(nghA[i])!=setB.end()) {
                ans += 1;
            }
        }
    }

    return ans;
}

template <class vertex>
struct countF{
    vertex* V;
    long* counts;
    countF(vertex*_V, long* _counts) : V(_V), counts(_counts) {}
    inline bool update(uintE s, uintE d) {
        if (s > d) {
            writeAdd(&counts[s], countCommon<vertex>(V[s],V[d],s,d));
        }
        return 1;
    }
    inline bool updateAtomic(uintE s, uintE d) {
        if (s>d) {
            writeAdd(&counts[s], countCommon<vertex>(V[s],V[d],s,d));
        }
        return 1;
    }
    inline bool cond(uintE d) {return cond_true(d);}
};

struct intLT {bool operator() (uintT a, uintT b){return a<b;}};

template<class vertex> 
struct initF {
    vertex* V;
    long * counts;
    initF(vertex* _V, long* _counts): V(_V), counts(_counts) {}
    inline bool operator () (uintE i) {
        counts[i] = 0;
        return 1;
    }
};

template <class vertex>
void Compute(graph<vertex> & GA, commandLine P) {
    uintT n = GA.n;
    long* counts = newA(long, n);
    bool * frontier = newA(bool, n);
    {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
    vertexSubset Frontier(n, n, frontier);
    vertexMap(Frontier, initF<vertex>(GA.getvertex(), counts));
    edgeMap(&GA, Frontier, countF<vertex>(GA.getvertex(),counts), -1, no_output);
    long count = sequence::plusReduce(counts, n);
    cout << "triangle count = " << count << endl;
    Frontier.del(); free(counts);
}