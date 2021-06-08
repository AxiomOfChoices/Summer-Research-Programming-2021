# WeightManifold

## A small class to extract the quotient manifold of weight conformal equivalence classes of an arbitrary graph

To use put in a folder of a sagemath file then use 
```
    from WeightManifold import WeightManifold
```
Then to get the manifold of a graph ```G``` use
```
    Manifold = WeightManifold(G)
```

## Documentation

```
    Manifold.get_dimension()
```
Returns integer with the dimension of the manifold

```
    Manifold.get_basis()
```
Returns matrix with basis vectors as column vectors

```
    Manifold.get_arbitrary_element()
```
Returns coordinate form of the unique weight vector in the conformal class [w] with the property that the product of edges around a vertex is always 1.

```
    Manifold.submanifold_equation()
```
Returns defining function of the submanifold with determinant of adjecency matrix equal to 0 