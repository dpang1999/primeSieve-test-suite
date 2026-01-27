import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';

export class DoubleField implements IField<DoubleField>, IMath<DoubleField>, IOrdered<DoubleField>{
  value: number;
  constructor(value: number) {
    this.value = value;
  }
  a(o: DoubleField): DoubleField { return new DoubleField(this.value + o.value); }
  ae(o: DoubleField): void { this.value += o.value; }
  s(o: DoubleField): DoubleField { return new DoubleField(this.value - o.value); }
  se(o: DoubleField): void { this.value -= o.value; }
  m(o: DoubleField): DoubleField { return new DoubleField(this.value * o.value); }
  me(o: DoubleField): void { this.value *= o.value; }
  d(o: DoubleField): DoubleField { return new DoubleField(this.value / o.value); }
  de(o: DoubleField): void { this.value /= o.value; }
  coerce(o: number): DoubleField { return new DoubleField(o); }
  coerce_to_number(): number { return this.value; }
  is_zero(): boolean { return this.value === 0; }
  is_one(): boolean { return this.value === 1; }
  zero(): DoubleField { return new DoubleField(0); }
  one(): DoubleField { return new DoubleField(1); }
  abs(): number { return Math.abs(this.value); }
  sqrt(): DoubleField { return new DoubleField(Math.sqrt(this.value)); }
  lt(o: DoubleField): boolean { return this.value < o.value; }
  le(o: DoubleField): boolean { return this.value <= o.value; }
  gt(o: DoubleField): boolean { return this.value > o.value; }
  ge(o: DoubleField): boolean { return this.value >= o.value; }
  eq(o: DoubleField): boolean { return this.value === o.value; }
  copy(): DoubleField { return new DoubleField(this.value); }
  toString(): string { return this.value.toFixed(4); }
}
