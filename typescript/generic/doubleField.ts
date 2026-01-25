import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { ICopiable } from './iCopiable';

export class DoubleField implements IField<DoubleField>, IMath<DoubleField>, IOrdered<DoubleField>, ICopiable<DoubleField> {
  value: number;
  constructor(value: number) {
    this.value = value;
  }
  add(o: DoubleField): DoubleField { return new DoubleField(this.value + o.value); }
  sub(o: DoubleField): DoubleField { return new DoubleField(this.value - o.value); }
  mul(o: DoubleField): DoubleField { return new DoubleField(this.value * o.value); }
  div(o: DoubleField): DoubleField { return new DoubleField(this.value / o.value); }
  coerceFromInt(i: number): DoubleField { return new DoubleField(i); }
  coerceFromFloat(f: number): DoubleField { return new DoubleField(f); }
  coerceToFloat(): number { return this.value; }
  isZero(): boolean { return this.value === 0; }
  isOne(): boolean { return this.value === 1; }
  zero(): DoubleField { return new DoubleField(0); }
  one(): DoubleField { return new DoubleField(1); }
  abs(): DoubleField { return new DoubleField(Math.abs(this.value)); }
  sqrt(): DoubleField { return new DoubleField(Math.sqrt(this.value)); }
  lt(o: DoubleField): boolean { return this.value < o.value; }
  le(o: DoubleField): boolean { return this.value <= o.value; }
  gt(o: DoubleField): boolean { return this.value > o.value; }
  ge(o: DoubleField): boolean { return this.value >= o.value; }
  eq(o: DoubleField): boolean { return this.value === o.value; }
  copy(): DoubleField { return new DoubleField(this.value); }
  toString(): string { return this.value.toString(); }
}
