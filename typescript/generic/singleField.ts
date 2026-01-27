import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';

export class SingleField implements IField<SingleField>, IMath<SingleField>, IOrdered<SingleField>{
  value: number;
  constructor(value: number) {
    this.value = value;
  }
  a(o: SingleField): SingleField { return new SingleField(this.value + o.value); }
  ae(o: SingleField): void { this.value += o.value; }
  s(o: SingleField): SingleField { return new SingleField(this.value - o.value); }
  se(o: SingleField): void { this.value -= o.value; }
  m(o: SingleField): SingleField { return new SingleField(this.value * o.value); }
  me(o: SingleField): void { this.value *= o.value; }
  d(o: SingleField): SingleField { return new SingleField(this.value / o.value); }
  de(o: SingleField): void { this.value /= o.value; }
  coerce(o: number): SingleField { return new SingleField(o); }
  coerce_to_number(): number { return this.value; }
  is_zero(): boolean { return this.value === 0; }
  is_one(): boolean { return this.value === 1; }
  zero(): SingleField { return new SingleField(0); }
  one(): SingleField { return new SingleField(1); }
  abs(): number { return Math.abs(this.value); }
  sqrt(): SingleField { return new SingleField(Math.sqrt(this.value)); }
  lt(o: SingleField): boolean { return this.value < o.value; }
  le(o: SingleField): boolean { return this.value <= o.value; }
  gt(o: SingleField): boolean { return this.value > o.value; }
  ge(o: SingleField): boolean { return this.value >= o.value; }
  eq(o: SingleField): boolean { return this.value === o.value; }
  copy(): SingleField { return new SingleField(this.value); }
  toString(): string { return this.value.toFixed(4); }
}
