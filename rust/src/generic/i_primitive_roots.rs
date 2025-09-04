pub trait IPrimitiveRoots<N> {
    fn primitive_root(&self, p: u128) -> N;
    fn pow(&self, exp: u128) -> N;
    fn precomputeRootsOfUnity(&self, n: u32, direction: i32) -> Vec<N>;
}