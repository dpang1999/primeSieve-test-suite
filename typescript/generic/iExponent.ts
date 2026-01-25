export interface IExponent<T> {
  add(o: T): T;
  sub(o: T): T;
  lcm(o: T): T;
  deg(): number;
  lexCompare(o: T): number;
  canReduce(o: T): boolean;
  equals(o: T): boolean;
}