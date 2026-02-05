import { IField } from './iField';
import { IMath } from './iMath';
import { LCG } from '../helpers/lcg';
import { IOrdered } from './iOrdered';

export function monteCarloPi<C extends IField<C> & IMath<C> & IOrdered<C>>(field: C, n: number): number {
  const rand = new LCG(12345, 1345, 16645, 1013904);
  let underCurve = 0;
  for (let i = 0; i < n; i++) {
    const x = field.coerce(rand.nextDouble());
    const y = field.coerce(rand.nextDouble());
    // x^2 + y^2 <= 1
    const x2 = x.m(x);
    const y2 = y.m(y);
    const sum = x2.a(y2);
    if (sum.le(field.one())) {
      underCurve++;
    }
  }
  return (underCurve / n) * 4.0;
}


if (require.main === module) {
  // Usage: tsx genMonteCarlo.ts [numSamples] [mode]
  // mode: 1 = SingleField, 2 = DoubleField, else IntModP
  const args = process.argv.slice(2);
  let numSamples = 1_000_000;
  let mode = 2;
  if (args.length > 0) numSamples = parseInt(args[0], 10) || numSamples;
  if (args.length > 1) mode = parseInt(args[1], 10) || mode;

  let pi: number;
  if (mode === 1) {
    // SingleField (float32)
    const { SingleField } = require('./singleField');
    const field = new SingleField(0.0);
    pi = monteCarloPi(field, numSamples);
  } else if (mode === 2) {
    // DoubleField (float64)
    const { DoubleField } = require('./doubleField');
    const field = new DoubleField(0.0);
    pi = monteCarloPi(field, numSamples);
  } else {
    // IntModP (modular arithmetic, not meaningful for pi, but for completeness)
    const { IntModP } = require('./intModP');
    IntModP.setModulus(1_000_000_007);
    const field = new IntModP(0);
    pi = monteCarloPi(field, numSamples);
  }
  console.log(`Pi is approximately: ${pi}`);
  console.log(`Num samples: ${numSamples}`);
  console.log(`RMS Error: ${Math.abs(Math.PI - pi)}`);
}
