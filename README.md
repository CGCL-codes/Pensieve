Pensieve: Skewness-Aware Version Switching for Efficient Graph Processing
======================

Pensieve is a skewness-aware multi-version graph processing system. Two factors contribute to the efficiency of Pensieve. First, Pensieve leverages a differentiated graph storage strategy that stores low degree vertices using copy-based scheme while stores high degree ones using delta-based scheme. Such a design achieves a good trade-off between storage cost and version switching time for multi-version graph processing. Second, the Pensieve graph storage exploits the time locality of graph version access and designs a novel last-root version switching scheme, which significantly improves the access efficiency for recent versions. We implement Pensieve on top of Ligra, and conduct comprehensive experiments to evaluate the performance of this design using large-scale datasets collected from real world systems. The results show that Pensieve substantially outperforms state-of-the-art designs in terms of memory consumption and version switching time.

Introduction
--------

Graphs are widely used to capture relations among entities. For example, social graphs represent interconnections and interactions among users, and Web graphs denote links between web pages. Since many practical issues are equivalent to graph problems, graph processing has attracted much research efforts. 

As a graph commonly evolves constantly (e.g., social links change in Facebook), it produces a large number of versions as time flies. A graph version corresponds to a snapshot at a specific time of the evolving graph. 
A multi-version graph processing system executes tasks on different versions of diverse graphs. For example, a single analyst may need to perform binary search in graph history and run algorithms to find a mutation. Therefore, a graph processing system should support arbitrary graph version switching.

Since historical graphs are often too large to fit into memory, traditional graph processing systems store different graph versions in persistent storage. They switch graph versions by simply loading the desired graph version from persistent storage into memory and replacing the previously processed version. 
The straightforward graph loading is costly for multi-version graph processing. Existing schemes to solve this problem fall into two types: copy-based schemes and delta-based schemes.

A copy-based scheme saves both the old and new data when a new version is generated. For example, Version Traveler([USENIX ATC 2016](https://www.usenix.org/conference/atc16/technical-sessions)). The system generates and stores the copies of two new versions of the graph. If all the versions are stored in memory, a graph processing system can switch among versions very fast since nearly no extra processing is needed except simple pointer redirecting. However, such a scheme is memory-inefficient due to high redundancy in storage. 

A delta-based scheme stores a basic full graph version and the delta logs of added/deleted edges, like ASGraph. It computes a target version by adding the relevant deltas into the basic version. Such a scheme is more memory-efficient than the copy-based scheme by avoiding redundant storage. However, a delta-based scheme introduces longer latency for version switching due to extra delta computation. Such a delta-based scheme achieves memory efficient multi-version graph storage at the cost of extra computation during version switching.

Our insight is that the presence of skewness in graphs greatly influences the system efficiency. 
First, the degrees of vertices are highly skewed. A small fraction of high degree vertices dominates the storage overhead of the copy-based scheme.
Second, the access frequencies of different graph versions are skew. Users are commonly more interested in recent content. Graphs are used to repersent all kinds of data.

we present Pensieve, a novel skewness-aware multi-version graph processing system. Two factors contribute to the efficiency of Pensieve. First, Pensieve provides differentiated multi-version graph storage organization. For high-degree vertices, Pensieve uses delta-based scheme to save a large amount of memory cost with slight computation. For low-degree vertices, Pensieve leverages copy-based scheme to accelerate version switching with affordable memory overhead. Second, by exploiting the access skewness, Pensieve stores the most recent version as the root version and performs version switching in a backward manner. The evolution of a graph commonly incurs much more edge additions than deletions. Therefore backward version switching mainly involves edge deletion, which previous system designs try to avoid because it either takes time to search deleted item or memory to store redundantly. To address this problem, Pensieve proposes a novel deletion-efficient delta computation scheme.

Structure of Pensieve
--------
![Pensieve Architecture](https://github.com/Pensieve-code/Pensieve/raw/master/tutorial/Pensieve_arch.png)

Pensieve provides multi-version graph storage for a multi-version graph processing system. 
Pensieve stores graph data and delta separately. Both graph storage and delta storage are divided into two part: high degree vertices in graph storage correspond to delta-based delta format and low degree vertices to copy-based delta. 
Pensieve has two control components, including vertex splitter and version controller. The vertex splitter processes new deltas and stores them in different delta storage, while the version controller captures the relationship among versions. When the graph processing system needs a specific version, the version controller responds to the request and dispatches the task to delta storage and graph storage to perform version switching.


Compilation and Run
--------

Compilation is done from within the apps/ directory. 

Compilers

* g++ &gt;= 5.3.0 

After the appropriate environment variables are set, to compile,
simply run

```
$ make 
```

The following commands cleans the directory:
```
$ make clean #removes all executables
$ make cleansrc #removes all executables and linked files from the ligra/ directory
```

To run a simple tast (like BFS):
```
./BFS -n 100 -s /path/to/graph/file
```

Before you run the code, please make sure you have the necessary delta files. Path of these files should be configured in main function of ligra/ligra.h

There are plenty of optional in the main function in ligra/ligra.h. Please activate them as needed.


Input Format for Pensieve
-----------

Pensieve currently support unweighted graphs.

The input format of unweighted graphs should be the following formats.

The adjacency graph format from the Problem Based Benchmark Suite
(http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html). The adjacency
graph format starts with a sequence of offsets one for each vertex,
followed by a sequence of directed edges ordered by their source
vertex. The offset for a vertex i refers to the location of the start
of a contiguous block of out edges for vertex i in the sequence of
edges. The block continues until the offset of the next vertex, or
the end if i is the last vertex. All vertices and offsets are 0 based
and represented in decimal. The specific format is as follows:

### Graph

AdjacencyGraph  
\<n\>  
\<m\>  
\<o0\>  
\<o1\>  
...  
<o(n-1)\>  
&lt;e0>  
&lt;e1>  
...  
&lt;e(m-1)>  

### Delta

DELTA_FILE  
<version_start> <version_end>  
\<add_number\>  
\<del_number\>  
<add_src0> <add_dst0>    
<add_src1> <add_dst1>  
...  
<add_srcx> <add_dstx>  
<del_src0> <del_dst0>      
<del_src1> <del_dst1>  
...  
<del_srcy> <del_dsty>

All file are represented as plain text.


Graph Applications
---------
Implementation files are provided in the apps/ directory: **BFS.C**
(breadth-first search), **BFS-Bitvector.C** (breadth-first search with
a bitvector to mark visited vertices), **BC.C** (betweenness
centrality), **Radii.C** (graph eccentricity estimation),
**Components.C** (connected components), **BellmanFord.C**
(Bellman-Ford shortest paths), **PageRank.C**, **PageRankDelta.C**,
**BFSCC.C** (connected components based on BFS), **MIS.C** (maximal
independent set), **KCore.C** (K-core decomposition), **Triangle.C**
(triangle counting), and **CF.C** (collaborative filtering).


Publications  
-------- 
If you want to know more detailed information, please refer to this paper:  
Tangwei Ying, Hanhua Chen, Hai Jin. [Pensieve: Skewness-Aware Version Switching for Efficient Graph Processing](https://dl.acm.org/doi/10.1145/3318464.3380590) in Proceedings of the 2020 International Conference on Management of Data (SIGMOD 2020), online conference [Portland, OR, USA], June 14-19, 2020.  
([bibtex](https://github.com/CGCL-codes/Pensieve/blob/master/Pensieve.bib))


Authors and Copyright
--------

Pensieve is developed in National Engineering Research Center for Big Data Technology and System, Cluster and Grid Computing Lab, Services Computing Technology and System Lab, School of Computer Science and Technology, Huazhong University of Science and Technology, Wuhan, China by Tangwei Ying(ytw@hust.edu.cn), Hanhua Chen (chen@hust.edu.cn) and Hai Jin (hjin@hust.edu.cn).

Copyright (C) 2019, [STCS & CGCL](grid.hust.edu.cn) and [Huazhong University of Science and Technology](www.hust.edu.cn).