from __future__ import annotations
from typing import List
from sage.graphs.graph import Graph
from sage.calculus.var import var
from sage.functions.log import exp, ln
from sage.modules.free_module_element import vector
from sage.rings.real_mpfr import RR
from sage.rings.integer_ring import ZZ


class ConformalClass:
    def __init__(self, parent: WeightManifold, g, weights: (List[float] | None) = None) -> None:
        if (type(g) != Graph):
            raise TypeError("g must be a graph")
        if (g.size() != len(weights)):
            raise ValueError("g must have as many edges as inputted")
        self.graph = g
        self.weights = weights if weights else g.edge_labels()
        self.parent = parent
        self.canonical_representative = self.parent.projection_matrix * vector(list(map(ln,self.weights)))
        self.space = ((ZZ**(self.graph.order())).hom(self.parent.incidence_matrix,(ZZ**(self.graph.size())))).image()
        self.dimension = self.space.dimension()
        self.basis = self.space.basis_matrix().T
        self.coordinates = var('z', n=self.dimension, latex_name='z', domain='real')
        
    def get_parent(self):
        return self.parent

    def get_arbitrary_element(self):
        return list(map(exp, self.canonical_representative+self.basis*vector(self.coordinates)))

    def get_arbitrary_element_graph(self):
        graph = self.graph.copy()
        graph.weighted(True)
        weights = self.get_arbitrary_element()
        edges = list(graph.edge_iterator())
        for edge, weight in zip(edges, weights):
            graph.set_edge_label(edge[0], edge[1], weight);
        return graph

    def get_random_element(self):
        return list(map(exp, self.canonical_representative+(self.basis.change_ring(RR)*(self.space.change_ring(RR).random_element()))))

    def get_random_element_graph(self):
        graph = self.graph.copy()
        graph.weighted(True)
        weights = self.get_random_element()
        edges = list(graph.edge_iterator())
        for edge, weight in zip(edges, weights):
            graph.set_edge_label(edge[0], edge[1], weight);
        return graph

class WeightManifold:
    def __init__(self, g):
        if (type(g) != Graph):
            raise TypeError("g must be a graph")
        graph = g.copy()
        graph.weighted(True)
        # We give all edges weights so that the weighted matrices are well defined
        edges = list(graph.edge_iterator())
        self.vars = var('y', n=len(edges), latex_name='y', domain='real');
        for edge, variable in zip(edges, self.vars):
            graph.set_edge_label(edge[0], edge[1], variable);
        self.graph = graph
        self.adjacency_matrix = graph.weighted_adjacency_matrix()
        self.incidence_matrix = graph.incidence_matrix(oriented=False)
        # From previous papers we know that the quotient manifold is isomorphic
        # to the kernel of the incidence_matrix
        self.quotient_manifold = self.incidence_matrix.right_kernel()
        self.dimension = self.quotient_manifold.dimension()
        self.basis = self.quotient_manifold.basis_matrix().T
        self.projection_matrix = (
                                    self.basis
                                    * (self.basis.T * self.basis).inverse()
                                    * self.basis.T
                                 )
        if (self.dimension > 0):
            self.coordinates = var('x', n=self.dimension, latex_name='x', domain='real')
        else:
            self.coordinates = []
        

    def get_canonical_arbitrary_representative(self):
        # Multiply basis by coordinate vector and then exponentiate
        # it to get actual weights
        return list(map(exp, self.basis*vector(self.coordinates)))

    def get_arbitrary_weight_class(self):
        representative = self.get_canonical_arbitrary_representative()
        return ConformalClass(self, self.graph, representative)

    def get_random_canonical_weight_representative(self):
        return list(map(exp, self.quotient_manifold.change_ring(RR).random_element()))

    def get_random_canonical_representative_graph(self):
        weights = list(map(exp, self.quotient_manifold.random_element()))
        graph = self.graph.copy()
        graph.weighted(True)
        edges = list(graph.edge_iterator())
        for edge, weight in zip(edges, weights):
            graph.set_edge_label(edge[0], edge[1], weight);
        return graph

    def get_random_weight_class(self):
        representative = self.get_random_canonical_weight_representative()
        return ConformalClass(self, self.graph, representative)
        

    def submanifold_equation(self):
        # We construct a substitution list that we apply to the determinant of the
        # Adjacency matrix to get it in terms of coordinates
        substitutions = []
        for name, value in zip(self.vars, self.get_canonical_arbitrary_representative()):
            substitutions.append(name==value)
        deter = self.adjacency_matrix.det().subs(substitutions)
        # We need to remember that our coordinates are the logs of the weights
        # and so we substitute the log in order to get the equation in terms
        # of actual weights 
        substitutions = []
        for name in self.coordinates:
            substitutions.append(name==ln(name))
        deter = deter.subs(substitutions)
        deter = deter.simplify_full()
        return (deter == 0)

    


