package generic;

interface IPrimitiveRoots<T extends IField<T>> {
    T primitiveRoot(int n);
    T pow(int exp);
    T[] precomputeRootsOfUnity(int n, int direction);
}
