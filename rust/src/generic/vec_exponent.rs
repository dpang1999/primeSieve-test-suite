use std::fmt;
use crate::generic::i_exponent::IExponent;
use crate::generic::i_ordered::IOrdered;
use std::hash::{Hash, Hasher};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct VecExponent{
    pub exponents:  Vec<u32>
}

impl VecExponent {
    pub fn new(exponents: Vec<u32>) -> Self {
        VecExponent { exponents }
    }
}

impl Hash for VecExponent {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for exp in &self.exponents {
            exp.hash(state);
        }
    }
}

impl IExponent for VecExponent {
    fn add(&self, o: &Self) -> Self {
        let result: Vec<u32> = self.exponents.iter().zip(&o.exponents).map(|(a, b)| a + b).collect();
        VecExponent{ exponents: result }
    }

    fn sub(&self, o: &Self) -> Self {
        let result: Vec<u32> = self.exponents.iter().zip(&o.exponents).map(|(a, b)| a - b).collect();
        VecExponent{ exponents: result }
    }

    fn lcm(&self, other: &Self) -> Self {
        let result: Vec<u32> = self
            .exponents
            .iter()
            .zip(&other.exponents)
            .map(|(a, b)| std::cmp::max(*a, *b))
            .collect();
        VecExponent{ exponents: result }
    }

    fn degree(&self) -> u32 {
        self.exponents.iter().sum()
    }

    fn lex_compare(&self, other: &Self) -> std::cmp::Ordering {
        for (a, b) in self.exponents.iter().zip(&other.exponents) {
            if a < b {
                return std::cmp::Ordering::Less;
            } else if a > b {
                return std::cmp::Ordering::Greater;
            }
        }
        std::cmp::Ordering::Equal
    }

    fn can_reduce(&self, divisor: &Self) -> bool {
        for (a, b) in self.exponents.iter().zip(&divisor.exponents) {
            if a < b {
                return false;
            }
        }
        true
    }
}


impl fmt::Display for VecExponent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let degree: u32 = self.degree();
        write!(f, "Degree: {:04X}, Exponents: ", degree)?;
        for exp in &self.exponents {
            write!(f, "{:02X} ", exp)?;
        }
        Ok(())
    }
}