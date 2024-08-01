# The Ansatz

After running `maxcut_qaoa.jl`,
    the `ansatz` contains everything you need to create a circuit in qiskit or whatever.

In particular, there are four relevant attributes:
- `ansatz.observable` - a "ScaledPauliVector" $H$
- `ansatz.generators` - a list of "ScaledPauliVector" generators $G_l$, one for each layer
- `ansatz.β_parameters` - a list of floats $β_l$, one for each layer
- `ansatz.γ_parameters` - a list of floats $γ_l$, one for each layer

The circuit, in math form, looks like:

$$U = \exp(-i β_L G_L) \exp(-i γ_L H) ... \exp(-i β_1 G_1) \exp(-i γ_1 H)$$

NOTE: when using the `ADAPT.Ansatz` sort of ansatz, the only relevant attributes are `ansatz.generators` and `ansatz.parameters`.

# The ScaledPauliVector Type

The "ScaledPauliVector" is, quite simply, a list of "ScaledPauli", representing a sum of Paulis.
Documentation for this object is in [PauliOperators.jl](https://github.com/nmayhall-vt/PauliOperators.jl).

In `ansatz.generators`, each list (say the one indexed by `ansatz.generators[i]`) has only a single element
    (so long as you are using the `qaoa_double_pool`).
Fetch it with `only(ansatz.generators[i])`.

In `ansatz.observable`, the list has more than one "ScaledPauli", but they all commute with one another,
    so to do the exponentiation, you _could_ just exponentiate each individual "ScaledPauli" in sequence,
    all scaled by the same $γ_l$.
I am certain there are better ways of doing this in hardware, but I leave that to you. :)

So for now, let's say we're just trying to exponentiate a single "ScaledPauli".
Let's call it `sc_pauli`.
It has two attributes:
- `coeff` - a float number $c$
- `pauli` - a "FixedPhasePauli" $P$

In `ansatz.generators` (so long as you are using the `qaoa_double_pool`),
    all the coeffients are actually 1.0 and you can ignore them,
    although they may not look it at first glance.
More on that in the next section.

# The FixedPhasePauli Type

Let's call our "FixedPhasePauli" `pauli`.
It has two attributes, `pauli.x` and `pauli.z`, both of which are integers.

The relevant information in these integers is their bitstring representations.
E.g.
```
bitstring(pauli.x)
```
will give a string like `0000000000000...0010010`.
Only the `N` least significant bits matter, where `N` is the number of qubits you're working with.

If there is a 1 in the bitstring representation in the q-th least signifcant bit,
    it means this Pauli word acts like a Pauli-X acting on the q-th qubit.

Meanwhile, if there is a 1 in the same bit of `pauli.z`,
    it means this Pauli word acts like a Pauli-Z acting on the q-th qubit.

It may well act like both: X first, then Z, gives $ZX = iY$.
Note, therefore, that a "FixedPhasePauli" is _not_ necessarily a Pauli with a phase of _one_.
The fixed phase is actually $i^k$, where $k$ is the number of bits which have a 1 in both `pauli.x` and `pauli.z`.
You will find that the `sc_pauli.coeff` in each "ScaledPauli" found within `ansatz.generators`
    simply serves to cancel out this fixed phase,
    which is why I said you can ignore it.
But also please note that `ansatz.observable` _does_ have some scaling going on,
    so don't ignore it there. :)

There is a very good chance that whatever software you are trying to use
    has a way of constructing Paulis directly from strings,
    which is probably much easier than trying to parse the tableau representation.
So, perhaps all I needed to show in this section was the following:
```
string(pauli)
```
That'll give the `N`-length string of "I", "X", "y", and "Z",
    which you can probably plug directly into qiskit or whatever.
(But, FYI, qiskit probably also uses the tableau representation under the hood!)

Last note, to beware of: that string used a lower-case "y" rather than the standard "Y".
That's because it is really representing the product $ZX=iY$, and we didn't want to confuse.
This distinction doesn't matter for `ansatz.generators`,
    since the `sc_pauli.coeff` is indeed canceling out the phase as discussed above.
It also doesn't matter for `ansatz.observable`,
    simply because the observable in maxcut QAOA is always diagonal - there should only be "I" and "Z".

# Summary
This line of code is a pretty good cheat-sheet for everything described above:
```
string(only(ansatz.generators[i]).pauli)
```