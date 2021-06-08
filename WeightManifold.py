from sage.graphs.graph import Graph
from sage.calculus.var import var
from sage.functions.log import exp, ln
from sage.modules.free_module_element import vector

class WeightsManifold:
    def __init__(self, g):
        if (type(g) != Graph):
            raise Error("g must be a graph")
        graph = g.copy()
        graph.weighted(True)
        # We give all edges weights so that the weighted matrices are well defined
        edges = list(graph.edge_iterator())
        self.Vars = var('y', n=len(edges), latex_name='x', domain='real');
        for edge, variable in zip(edges, self.Vars):
            graph.set_edge_label(edge[0], edge[1], variable);
        self.graph = graph
        self.Adjacency_matrix = graph.weighted_adjacency_matrix()
        self.Incidence_matrix = graph.incidence_matrix()
        # From previous papers we know that the quotient manifold is isomorphic
        # to the kernal of the incidence_matrix
        self.QuotientManifold = self.Incidence_matrix.right_kernel()
        if (self.get_dimension() > 0):
            self.coordinates = var('x', n=self.get_dimension(), latex_name='x', domain='real')
        else:
            self.coordinates = []
        
    def get_dimension(self):
        return self.QuotientManifold.dimension()
    
    def get_basis(self):
        # Transpose so that we can multiply by 
        # coordinate vector to get arbitrary element
        return self.QuotientManifold.basis_matrix().transpose()
    
    def get_arbitrary_element(self):
        # Multiply basis by coordinate vector and then exponentiate
        # it to get actual weights
        return list(map(exp, self.get_basis()*vector(self.coordinates)))
    
    def submanifold_equation(self):
        # We construct a substitution list that we apply to the determinant of the
        # Adjacency matrix to get it in terms of coordinates
        substitutions = []
        for name,value in zip(self.Vars, self.get_arbitrary_element()):
            substitutions.append(name==value)
        deter = self.Adjacency_matrix.det().subs(substitutions)
        # We need to remember that our coordinates are the logs of the weights
        # and so we substitute the log in order to get the equation in terms
        # of actual weights 
        substitutions = []
        for name in self.coordinates:
            substitutions.append(name==ln(name))
        deter = deter.subs(substitutions)
        deter = deter.simplify_full()
        return (deter == 0)