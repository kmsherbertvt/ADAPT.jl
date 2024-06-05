# ADAPT

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kmsherbertvt.github.io/ADAPT.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kmsherbertvt.github.io/ADAPT.jl/dev/)
[![Build Status](https://github.com/kmsherbertvt/ADAPT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kmsherbertvt/ADAPT.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kmsherbertvt/ADAPT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kmsherbertvt/ADAPT.jl)

This is my implementation of [ADAPT-VQE](https://www.nature.com/articles/s41467-019-10988-2).
There are many like it, but this one is mine.

#### Why This Implementation

I like to write code. That's the main reason this one exists. ^_^

In particular, it is meant to address three things I found unsatisfactory in the existing Julia implementation:
1. More thorough awareness of the VQE optimization, ie. whether it "converges" and so on.
2. A "history" that can include a trace over all the optimizations, not just the adaptations.
3. A mechanism for, if it's ever necessary, working with generators that are not necessarily a sum of *commuting* Pauli words.

It is designed to be completely modular and extensible,
    such that different variants of ADAPT can usually be implemented by dropping in a single file,
    plus some very light changes to the main package to include the contents of the file in the package.

#### Where is the Documentation?

At the top of the readme, there is a badge labeled "Dev".

If you click on it, you'll find online documentation that is kept up to date with continuous integration.

I don't know how to use GitHub well enough to know how to get something to happen when you click the "Stable" badge.
That that link is broken should _not_ be interpreted to mean that the code is in some way unstable. ^_^

Anyway, while my overall structure can sometimes get a touch convoluted (the price of modularity and extensibility),
    my docstrings are usually pretty thorough and my code itself is usually relatively readable,
    so it should be a net positive experience. ^_^

#### Overall Structure

Thorough tutorial documentation does not yet exist for this package,
    though looking through scripts in the `test` folder should be edifying.

For now, here's a very brief orientation:

The "main" function in `ADAPT.jl` is `ADAPT.run!()`. The [documentation](https://kmsherbertvt.github.io/ADAPT.jl/dev/#ADAPT.run!-Tuple%7BAbstractAnsatz,%20Dict%7BSymbol,%20Any%7D,%20AdaptProtocol,%20OptimizationProtocol,%20AbstractVector,%20Any,%20Any,%20AbstractVector%7B%3C:AbstractCallback%7D%7D) shows it takes the following arguments:
- `ansatz` - this represents the algorithm's ["state"](https://en.wikipedia.org/wiki/Finite-state_machine) - typcially a sequence of generators ADAPT has selected, together with their variational parameters.
- `trace` - this is just a `Dict` keeping track of useful quantities at each iteration.
- `adapt` - this represents the algorithm for selecting a new generator.
- `vqe` - this represents the algorithm for optimizing a given sequence of generators.
- `pool` - this is a `Vector` of all the possible generators ADAPT can choose from - typically a sequence of commuting Pauli operators.
- `observable` - this represents the algorithm for "measuring" a scalar (e.g. energy) from a quantum state - typically a weighted sum of Pauli operators.
- `reference` - this represents the initial ["state"](https://en.wikipedia.org/wiki/Quantum_state) of the quantum computer, prior to applying the ansatz.
- `callbacks` - this is a `Vector` of objects representing extra steps to do at each iteration - typical runs will want at minimum a `Tracer` and a `ScoreStopper`.

The code in `src/core` defines the, let's say, _schemae_ and _relationships_ for each of these arguments, but no concrete objects.

The code in `src/base` defines concrete objects for the most typical use cases.
Even for atypical use cases, _most_ of these can use the "basic" definitions,
    defined in precisely the way you'll find in the example scripts.

The fun part is figuring out the thing which is unique about your implementation,
    and deciding whether or not you need to implement a new type,
    and deciding which functions need to be overridden for your new type.