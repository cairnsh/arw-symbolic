from sympy import *
from defs import *
from transient import stabilize
from transient import transient_from_data

class graph:
    def __init__(self, n, edges, sink):
        edgelist = set()
        for (a,b) in edges:
            if a < b:
                edgelist.add((a, b))
            else:
                edgelist.add((b, a))
        edges = {a: [] for a in range(n)}
        for (a,b) in edgelist:
            edges[a].append(b)
            edges[b].append(a)
        
        for vv in sink:
            edges[vv] = []
        
        def nf():
            return list(range(n))
        def ef(v):
            return edges[v]
        def sf():
            return sink
        self.nodes = nf
        self.edges = ef
        self.sink = sf

class utils:
    @staticmethod
    def issink(position):
        return not any(z >= 1 for z in position)
    
    @staticmethod
    def moves(graph, sleeprates, position):
        if any(z > 1 for z in position):
            return subproblem.solve(graph, position)
        moves = utils.moves_unstabilized(graph, sleeprates, position)
        to = {}
        def add(pos, z):
            to[pos] = to.get(pos, 0) + z
        for move, prob in moves.items():
            print(move)
            stab = subproblem.solve(graph, move)
            for move2, prob2 in stab.items():
                add(move2, prob * prob2)
        return to
    
    @staticmethod
    def moves_unstabilized(graph, sleeprates, position):
        def moves_when_we_topple(j):
            to = graph.edges(j)
            n = len(to)
            if n == 0:
                return {position: 1}

            assert position[j] == 1
            sleep = position[:j] + (STOP,) + position[j+1:]
            moves = {sleep: sleeprates[j]}
            for i in to:
                pos = list(position)
                pos[j] -= 1
                if pos[i] == SINK:
                    pass
                elif pos[i] == STOP:
                    pos[i] = 2
                else:
                    pos[i] += 1
                pos = tuple(pos)
                if pos in moves:
                    raise Exception("shouldn't happen")
                moves[pos] = (1 - sleeprates[j]) * FRAC(1, n)
            return moves

        # We choose the unstable site with the lowest
        # index and topple it. (We could use a different
        # rule, although we'd get the same results in
        # the end by the abelian property)
        for j in range(len(position)):
            if position[j] != STOP and position[j] != 0 and position[j] != SINK:
                return moves_when_we_topple(j)
        return {position: 1}

class subproblem:
    # find the transition from a state with a 2 in it to states with only s and 1

    @staticmethod
    def sinkify(graph, a):
        foo = list(a)
        for v in graph.sink():
            foo[v] = SINK
        return tuple(foo)

    @staticmethod
    def edges(graph, a):
        a = subproblem.sinkify(graph, a)
        idx = None
        # find the two and make sure there's only one
        for j in range(len(a)):
            if a[j] == 2:
                if idx is None:
                    idx = j
                else:
                    raise Exception("more than one two in configuration %s" % str(a))
            elif a[j] > 2:
                raise Exception("more than two chips in configuration %s" % str(a))
        if idx is None:
            return []
        configs = []
        for other in graph.edges(idx):
            cur = list(a)
            cur[idx] = 1
            if cur[other] == STOP:
                cur[other] = 1
            if cur[other] == SINK:
                pass
            else:
                cur[other] += 1
            configs.append(tuple(cur))
        return configs
    
    @staticmethod
    def create_subproblem(graph, a):
        a = subproblem.sinkify(graph, a)
        nodes = []
        sink = []
        edges = {}
        frontier = [a]
        seen = {a: True}
        while frontier:
            what = frontier.pop()
            nodes.append(what)
            edges[what] = subproblem.edges(graph, what)
            if not edges[what]:
                sink.append(what)
            for vv in edges[what]:
                if vv in seen:
                    pass
                else:
                    frontier.append(vv)
                    seen[vv] = True
        def even(z):
            n = len(z)
            the = {}
            for j in z:
                the[j] = the.get(j, 0) + FRAC(1, n)
            return the
        edges = {a:even(b) for (a,b) in edges.items()}
        return nodes, edges, sink
    
    def solve(graph, a):
        the = subproblem.create_subproblem(graph, a)
        nodes, edges, sink = the
        if not sink:
            raise Exception("no stable configuration can be reached")
        print("stabilize", the, a)
        (out,) = stabilize(*the, [a])
        return out

def problem(graph, sleeprates, starting_position):
    # fill in the sink vertices with the sink indicator (-404)
    starting_position = subproblem.sinkify(graph, starting_position)
    boundary = [starting_position]
    visited = {}
    nodes = []
    edges = {}
    sink = []
    while boundary:
        position = boundary.pop()
        visited[position] = True
        nodes.append(position)
        edges[position] = utils.moves(graph, sleeprates, position)
        print("edges[%s] = %s" % (position, edges[position]))
        if utils.issink(position):
            sink.append(position)
        for node in edges[position]:
            if node not in visited:
                boundary.append(node)
    return nodes, edges, sink