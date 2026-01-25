export interface IPrimitiveRoots<T> {
  primitive_root(p: number): T;
  pow(exp: number): T;
  precomputeRootsOfUnity(n: number, direction: number): T[];
}
