import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { ICopiable } from './iCopiable';

export class SingleField implements IField<SingleField>, IMath<SingleField>, IOrdered<SingleField>, ICopiable<SingleField> {
  value: number;
  constructor(value: number) {
    this.value = value;
  }
  add(o: SingleField): SingleField { return new SingleField(this.value + o.value); }
  sub(o: SingleField): SingleField { return new SingleField(this.value - o.value); }
  mul(o: SingleField): SingleField { return new SingleField(this.value * o.value); }
  div(o: SingleField): SingleField { return new SingleField(this.value / o.value); }
  coerceFromInt(i: number): SingleField { return new SingleField(i); }
  coerceFromFloat(f: number): SingleField { return new SingleField(f); }
  coerceToFloat(): number { return this.value; }
  isZero(): boolean { return this.value === 0; }
  isOne(): boolean { return this.value === 1; }
  zero(): SingleField { return new SingleField(0); }
  one(): SingleField { return new SingleField(1); }
  abs(): SingleField { return new SingleField(Math.abs(this.value)); }
  sqrt(): SingleField { return new SingleField(Math.sqrt(this.value)); }
  lt(o: SingleField): boolean { return this.value < o.value; }
  le(o: SingleField): boolean { return this.value <= o.value; }
  gt(o: SingleField): boolean { return this.value > o.value; }
  ge(o: SingleField): boolean { return this.value >= o.value; }
  eq(o: SingleField): boolean { return this.value === o.value; }
  copy(): SingleField { return new SingleField(this.value); }
  toString(): string { return this.value.toString(); }
}
