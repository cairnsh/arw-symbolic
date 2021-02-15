computing arw steady state as a rational polynomial in the sleep rates

# Documentation

`import transient`

To make a transient:

```
class transient(abstract_transient):
    def nodes(self):
        # return node names
    
    def edges(self, vertex):
        # return edges from that vertex in
        # the form {vertex: probability}

    def sink(self):
        # return sink names

transient_from_graph(V, E, sink)
# makes a transient from a graph.
# V is list of names, E is list of pairs, and
# the nodes listed in sink are made absorbing.

stabilize(V, E, sink, startingnodes)
# return the stable distribution of the transient here
```
