pub trait IOrdered {
    fn lt (&self, o: &Self) -> bool;
    fn le (&self, o: &Self) -> bool;
    fn gt (&self, o: &Self) -> bool;
    fn ge (&self, o: &Self) -> bool;
    fn e (&self, o: &Self) -> bool;
}