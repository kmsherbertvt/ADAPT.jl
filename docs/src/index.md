```@meta
CurrentModule = ADAPT
```

# ADAPT

Documentation for [ADAPT](https://github.com/kmsherbertvt/ADAPT.jl).

```@index
```

### Core
```@autodocs
Modules = [
    ADAPT,
]
```

### Basics
```@autodocs
Modules = [
    ADAPT.Basics,
    ADAPT.Basics.Callbacks,
    ADAPT.Basics.Operators,
    ADAPT.Basics.Pools,
]
```

### Other Modules
```@autodocs
Modules = [
    ADAPT.OptimizationFreeADAPT,
    ADAPT.OverlapADAPT,
    ADAPT.Degenerate_ADAPT,
    ADAPT.Hamiltonians,
]
```

### MyPauliOperators
These methods should not be considered part of "ADAPT",
    but rather, destined for the `PauliOperators.jl` package.
The only reason I document them here is that
    the doc builder is configured to throw an error
    if any doc strings aren't included in the documentation...
```@autodocs
Modules = [
    ADAPT.Basics.MyPauliOperators,
]
```