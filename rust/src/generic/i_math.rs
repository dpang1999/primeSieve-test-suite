use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;
pub trait IMath: IField + IOrdered {
    fn abs(&self) -> Self;
    fn sqrt(&mut self);
}