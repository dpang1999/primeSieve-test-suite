use crate::generic::i_ring::IRing;
use crate::generic::i_ordered::IOrdered;
pub trait IMath: IRing + IOrdered {
    fn abs(&self) -> Self;
    fn sqrt(&mut self);
}