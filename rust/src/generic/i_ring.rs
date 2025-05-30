pub trait IRing: {
    fn a(&self, o: &Self) -> Self;
    fn ae(&mut self, o: &Self);
    fn s(&self, o: &Self) -> Self;
    fn se(&mut self, o: &Self);
    fn m(&self, o: &Self) -> Self;
    fn me(&mut self, o: &Self);

    fn coerce_from_i32(&self, i: i32) -> Self;
    fn coerce_from_f64(&self, i: f64) -> Self;
    fn coerce(&self) -> f64;

    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn zero() -> Self;
    fn one() -> Self;

    fn copy(&self) -> Self;
}