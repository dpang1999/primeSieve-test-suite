import { IField } from './iField';
import { IMath } from './iMath';

// Generic Monte Carlo integration for functions over any field
export function genericMonteCarlo<C>(
  f: (x: C) => C,
  a: C,
  b: C,
  field: IField<C> & IMath<C>,
  n: number
): C {
  let sum = field.zero();
  for (let i = 0; i < n; i++) {
    const t = Math.random();
    const x = field.add(a, field.mul(field.sub(b, a), field.fromNumber(t)));
    sum = field.add(sum, f(x));
  }
  return field.div(sum, field.fromNumber(n));
}
