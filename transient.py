# HAC 2020-10-15

from sympy import *
from progressbar import progressbar
import time
import os
from defs import FRAC
from solver import inv
import collections

class abstract_transient:
    """
    An interface class for a transient markov chain.
    See the code of print_transitions for an example.
    """
    def nodes(self):
        "Return a list of immutable node names."
        raise NotImplementedError()
    
    def edges(self, vertex):
        "return dictionary of transitions"
        raise NotImplementedError()
    
    def sink(self):
        "return list of sinks"
        raise NotImplementedError()

class transient_from_graph(abstract_transient):
    def __init__(self, V, E, sink):
        self.__nodes = V
        self.__sink = sink
        sinkset = set(sink)
        edgelist = {a: [] for a in V}
        for (a,b) in E:
            edgelist[a].append(b)
            edgelist[b].append(a)
        degree = {a: len(edgelist[a]) for a in V}
        def edges_for(a):
            if a in sinkset:
                return {}
            return {b: FRAC(1, degree[a]) for b in edgelist[a]}
        self.__edges = {a: edges_for(a) for a in V}
    
    nodes = lambda self: self.__nodes
    edges = lambda self, vertex: self.__edges[vertex]
    sink = lambda self: self.__sink

def print_transitions(tr):
    sink = tr.sink()
    sinkset = set(sink)
    for a in tr.nodes():
        edges = tr.edges(a)
        if len(edges) == 0:
            if a in sinkset:
                pass # we will print the sinks later
            else:
                raise Exception("node %s has no outgoing edges but is not a sink" % a)
        else:
            if a in sinkset:
                raise Exception("node %s has outgoing edges but is a sink" % a)
            print("node %s has transitions" % a)
            total = 0
            for (b,p) in edges.items():
                total += p
                print("    %-24s p = %s" % (b, p))
            assert simplify(total - 1) == 0
            print()
    print("sinks:")
    for a in sink:
        print("    %s" % a)

class reversible:
    def __init__(self, data=[]):
        self.data, self.reve = data, {data[i]: i for i in range(len(data))}
    def __getitem__(self, j):
        return self.data[j]
    def __len__(self):
        return len(self.data)
    def reverse(self, c):
        return self.reve[c]
    def has(self, c):
        return c in self.reve
    def __str__(self):
        return "reversible([%s])" % ", ".join([str(z) for z in self.data])
    def __repr__(self):
        out = ["%d: %s" % (i, repr(z)) for i, z in enumerate(self.data)]
        return "reversible([%s])" % ", ".join(out)
    
class transient_representation:
    "Represent a transient chain as a matrix like A_{ij} = P_i(X_1 = j)."
    def __init__(self, transient):
        self.nodes = reversible(transient.nodes())
        self.sink = reversible(transient.sink())
        self.transient = transient
        self.repr = 0
        
    def _get(self):
        if self.repr:
            return self.repr
        
        N = len(self.nodes)
        A = SparseMatrix(N, N, {})
        progress_bar = N>=100
        if progress_bar:
            pb = progressbar()
            t = time.time()

        entries = 0

        for j in range(N):
            if progress_bar and time.time() > t + 0.1:
                pb.display("line %d/%d" % (j, N))
                t = time.time()
            edges = self.transient.edges(self.nodes[j])
            total = sum(edges.values())
            if total != 0:
                for (b,p) in edges.items():
                    A[j, self.nodes.reverse(b)] += p/total
                    entries += 1
            else:
                assert self.sink.has(self.nodes[j])

        if progress_bar:
            pb.display("%d/%d complete" % (N, N))
            pb.done()
            print("constructed a sparse matrix with", entries, "entries")
            print("sparsity: %.4f%%" % (100 * entries / N / N))
            print("average entries per row: %.2f" % (entries / N))
        self.repr = A
        return A

    def get(self, coeff_of_identity, coeff_of_matrix):
        aa = self._get() * coeff_of_matrix
        for i in range(aa.rows):
            aa[i, i] += coeff_of_identity
        return aa
    
    def stable(self, startingnodes):
        b = SparseMatrix(len(startingnodes), len(self.nodes), {})
        for (i,c) in enumerate(startingnodes):
            b[i, self.nodes.reverse(c)] = 1
        # row i of b has 1 in the column self.node2ix[c] and 0 everywhere else
        a = self.get(1, -1)
        # xa = b
        x = inv(a, b)
        # extract the parts corresponding to sink
        sink_indices = [self.nodes[b] for b in self.sink]
        def stableof(i):
            return [x[i, j] for j in sink_indices]
        return [stableof(i) for i in range(len(startingnodes))]

def outerboundary(graph, region):
    outerboundary = []
    nodelist = set(region)
    for c in region:
        for e in graph.edges(c):
            if e not in nodelist:
                outerboundary.append(e)
                nodelist.add(e)
    return outerboundary

class linetransient(abstract_transient):
    def __init__(self, L, p):
        self.L = L
        self.p = p

    def nodes(self):
        return [j for j in range(self.L + 1)]
    
    def edges(self, j):
        if j == 0 or j == self.L: return {}
        return {j-1:1-self.p, j+1:self.p}

    def sink(self):
        return [0, self.L]
