export interface IOrdered<T> {
  lt(o: T): boolean;
  le(o: T): boolean;
  gt(o: T): boolean;
  ge(o: T): boolean;
  e(o: T): boolean;
}
