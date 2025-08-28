pub trait IPrimitiveRoots<N> {
    fn primitive_root(&self, p: u64) -> N;
    fn pow(&self, exp: i32) -> N;
    fn precomputeRootsOfUnity(&self, n: u64, direction: i32) -> Vec<N>;
}