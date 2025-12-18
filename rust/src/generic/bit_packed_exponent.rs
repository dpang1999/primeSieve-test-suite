use crate::generic::i_exponent::IExponent;
use std::hash::{Hash, Hasher};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BitPackedExponent {
    pub exponents: u64
}
impl BitPackedExponent {
    pub fn new(exponents: u64) -> Self {
        BitPackedExponent { exponents }
    }

    pub fn from_vec(exponents: [u8; 6]) -> Self {
        let mut packed: u64 = 0;
        
        // Pack the exponents with the first exponent at the left (bits 47..40)
        // Layout: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5
        for (i, &exp) in exponents.iter().enumerate() {
            let shift = 40_u32.saturating_sub((8 * i) as u32);
            packed |= (exp as u64) << shift;
        }

        // Compute the total degree (sum of all exponents)
        let degree: u16 = exponents.iter().map(|&e| e as u16).sum();

        // Add the degree to the top 16 bits of the u64
        packed |= (degree as u64) << 48;

        BitPackedExponent { exponents: packed }
    }
}

impl Hash for BitPackedExponent {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.exponents.hash(state);
    }
}

impl IExponent for BitPackedExponent {
    fn add(&self, o: &Self) -> Self {
        BitPackedExponent { exponents: self.exponents + o.exponents }
    }

    fn sub(&self, o: &Self) -> Self {
        BitPackedExponent { exponents: self.exponents - o.exponents }
    }

    fn lcm(&self, other: &Self) -> Self {
        let self_exponents = self.exponents & 0x0000_FFFF_FFFF_FFFF;
        let other_exponents = other.exponents & 0x0000_FFFF_FFFF_FFFF;

        let mut lcm_exponents: u64 = 0;
        let mut degree: u16 = 0;

        for i in (0..48).step_by(8) {
            let self_exp = (self_exponents >> i) & 0xFF;
            let other_exp = (other_exponents >> i) & 0xFF;

            let lcm_exp = self_exp.max(other_exp);

            lcm_exponents |= lcm_exp << i;
            degree += lcm_exp as u16; // Accumulate the degree
        }

        // Set the degree in the top 16 bits
        lcm_exponents |= (degree as u64) << 48;

        BitPackedExponent { exponents: lcm_exponents }
    }

    fn degree(&self) -> u32 {
        ((self.exponents >> 48) & 0xFFFF) as u32
    }

    fn lex_compare(&self, other: &Self) -> std::cmp::Ordering {
        for i in (0..48).step_by(8).rev() {
            let self_exp = self.exponents & 0x0000_FFFF_FFFF_FFFF;
            let other_exp = other.exponents & 0x0000_FFFF_FFFF_FFFF;

            if self_exp < other_exp {
                return std::cmp::Ordering::Less;
            } else if self_exp > other_exp {
                return std::cmp::Ordering::Greater;
            }
        }
        std::cmp::Ordering::Equal
    }

    fn can_reduce(&self, divisor: &Self) -> bool {
        for i in (0..48).step_by(8) {
            let self_exp = (self.exponents >> i) & 0xFF;
            let divisor_exp = (divisor.exponents >> i) & 0xFF;

            if self_exp < divisor_exp {
                return false;
            }
        }
        true
    }

}
