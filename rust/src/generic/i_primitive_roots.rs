pub trait IPrimitiveRoots<N> {
    fn primitive_root(&self, p: u64) -> N;
    fn pow(&self, exp: u64) -> N;
    fn precomputeRootsOfUnity(&self, n: u32, direction: i32) -> Vec<N>;
}