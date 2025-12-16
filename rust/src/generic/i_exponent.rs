pub trait IExponent {
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;
    fn degree(&self) -> u32;
    fn lex_compare(&self, other: &Self) -> std::cmp::Ordering;
    fn can_reduce(&self, divisor: &Self) -> bool;

}