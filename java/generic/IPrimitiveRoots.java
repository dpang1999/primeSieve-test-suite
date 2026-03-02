package generic;

interface IPrimitiveRoots<T> {
    T primitiveRoot(long n);
    T pow(long exp);
    T[] precomputeRootsOfUnity(int n, int direction);
}
