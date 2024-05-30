"""

This script won't actually run!!!

It's just an example stepping through how one might choose
    to plug in observables from `openfermion`,
    which can be particularly helpful when interfacing with,
    say, `pyscf` to generate molecular Hamiltonians.

We can state the problem thus:
    we have a Python object, a `QubitOperator`, in `openfermion`,
    which represents an observable we want to plug into our Julia ADAPT code.
The Julia ADAPT code expects, say, a `PauliSum` from the `PauliOperators.jl` package.
They represent the same information, but we need to translate from one to the other.


There's more than one way to interface with `openfermion`, of course.

For example, you could use `PyCall.jl` to work with `openfermion` data structures directly.
Then you could write a Julia function
    which takes the `QubitOperator` and returns a `PauliSum`.
But that requires having Julia set up within an appropriate Python environment,
    which I find tedious and undesirable.

My personal (convoluted) choice is to run Python separately
    to generate an `openfermion` data structure,
    and serialize it to a JSON file format.
Then I run my Julia code, load the JSON file, and parse it into a `PauliSum`.

This script has the code I've used to do that. (Not necessarily the most sensible!)

"""

the_python_code = """
    import json
    import openfermion

    def write_to_json(qubitop, filename):
        terms = [
            { "pauli": key, "coeff": value }
                for key, value in qubitop.terms.items()
        ]

        with open(f"{filename}.json", "w") as fp:
            json.dump(terms, fp)

    qubitop = # Do something with openfermion here!
    write_to_json(qubitop, "example")

"""

import PauliOperators
import PauliOperators: PauliSum, Pauli

""" Convert an json-form Pauli object into a dict-form Pauli object.

JSON-form:
[ {"pauli"=>[ [q, "X|Y|Z" ], "coeff"=>c ]

Dict-form:
{ (X=[q], Y=[q], Z=[q]) => c }

"""
function as_dict(json)
    dict = Dict{Any,Float}()

    for term in json
        # COLLECT ALL INDICES
        indices = Dict("X"=>Int[], "Y"=>Int[], "Z"=>Int[])
        for (q, op) in term["pauli"]
            push!(indices[op], 1+q) # Add one to convert from Python indexing to Julia
        end

        # SORT ALL INDICES
        # NOTE: pry redundant, if json came from openfermion, but req. for uniqueness
        foreach(op -> sort!(indices[op]), keys(indices))

        # CONVERT INDICES TO IMMUTABLE TUPLE
        key = (
            X = tuple(indices["X"]...),
            Y = tuple(indices["Y"]...),
            Z = tuple(indices["Z"]...),
        )

        dict[key] = term["coeff"]
    end

    return dict
end

""" Convert a Dict-form Pauli object into a PauliSum. """
function as_paulisum(n, dict)
    # CONSTRUCT A FUNCTION `qubitmap` SUCH THAT qubitmap(q) GIVES THAT QUBIT'S INDEX
    allqubits = Int[]
    for key in keys(dict)
        support = union(key.X, key.Y, key.Z)
        union!(allqubits, support)
    end
    sort!(allqubits)
    QUBITMAP = Dict(q => i for (i, q) in enumerate(allqubits))
    qubitmap = q -> QUBITMAP[q]

    # CONSTRUCT THE PAULI SUM
    pauli = PauliSum(n)
    for (key, coeff) in pairs(dict)
        mapped = (
            X = qubitmap.(key.X),
            Y = qubitmap.(key.Y),
            Z = qubitmap.(key.Z),
        )
        sum!(pauli, coeff * Pauli(n; mapped...))
    end
    return pauli
end


##########################################################################################
#= The actual script. =#

import JSON

#= Have already created a json file, using the code given in `the_python_code`. =#

# Load the JSON file, written in a form intuitively suited to `openfermion.QubitOperator`.
json = JSON.parsefile("example.json")

# Rearrange the data into a new dict inutiively suited to `PauliOperators.PauliSum`.
dict = MolecularADAPT.as_dict(json)

# Create the `PauliSum` object, which can be used as the "observable" in ADAPT.jl.
observable = as_paulisum(n, dict)