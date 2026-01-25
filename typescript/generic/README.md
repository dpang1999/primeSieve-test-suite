# TypeScript Generic Grobner Basis Implementation

This directory contains generic, modular, and idiomatic TypeScript implementations of Grobner basis algorithms, supporting multiple coefficient and exponent types. The design is inspired by the Rust and Go versions, with a focus on extensibility and testability.

## Key Files

- `field.ts`, `math.ts`, `ordered.ts`, `copied.ts`, `primitiveRoots.ts`: Core interfaces for fields, math, ordering, copying, and primitive roots.
- `doubleField.ts`, `singleField.ts`, `intModP.ts`: Field implementations for double, single, and modular arithmetic.
- `vecExponents.ts`, `bitPackedExponents.ts`: Exponent representations (vector and bitpacked).
- `polynomial.ts`: Generic `Term` and `Polynomial` classes.
- `grobner.ts`: Generic Grobner basis algorithm and polynomial operations.
- `genGrobner.ts`: Argument-driven random basis/test harness, supporting all combinations of coefficient and exponent types.
- `testGenGrobner.ts`: Test harness that runs all combinations of field and exponent types, and all term orders.

## Usage

To run the test harness (after compiling with a TypeScript compiler):

```
ts-node typescript/generic/testGenGrobner.ts
```

Or compile to JavaScript and run with Node.js:

```
tsc typescript/generic/testGenGrobner.ts --esModuleInterop
node typescript/generic/testGenGrobner.js
```

## Extending

- Add new field or exponent types by implementing the relevant interfaces and updating the random basis generator in `genGrobner.ts`.
- Specialized and additional algorithms (FFT, LU, SOR, Monte Carlo) can be added following this generic pattern.
