# ADAPT

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kmsherbertvt.github.io/ADAPT.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kmsherbertvt.github.io/ADAPT.jl/dev/)
[![Build Status](https://github.com/kmsherbertvt/ADAPT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kmsherbertvt/ADAPT.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kmsherbertvt/ADAPT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kmsherbertvt/ADAPT.jl)

This repo is an implementation of ADAPT-VQE
    meant to achieve three things I found unsatisfactory in the existing Julia implementation:
1. More thorough awareness of the VQE optimization, ie. whether it "converges" and so on.
2. A "history" that can include a trace over all the optimizations, not just the adaptations.
3. A mechanism for, if it's ever necessary, working with generators that are not necessarily a sum of *commuting* Pauli words.

It's also an excuse to use the Pauli library that Nick and Diksha and Griffin and Robert in the chemistry group have been tinkering with.

Ideally it will be sufficiently extensible to cover all the variants of ADAPT we've ever thought of.
But that's not the main goal, and it might not be possible, anyways... ^_^