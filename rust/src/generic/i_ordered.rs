use crate::generic::i_ring::IRing;
pub trait IOrdered: IRing {
    fn lt (&self, o: &Self) -> bool;
    fn le (&self, o: &Self) -> bool;
    fn gt (&self, o: &Self) -> bool;
    fn ge (&self, o: &Self) -> bool;
}