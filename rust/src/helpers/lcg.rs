pub struct Lcg {
    m: i32,
    a: i32,
    c: i32,
    last_num: i32,
}

impl Lcg {
    pub fn new(seed: i32, m: i32, a: i32, c: i32) -> Self {
        Self {
            m,
            // For the sake of fairness, modulus will always be 2^32 to bypass modulus bias between languages
            a,
            c,
            last_num: seed,
        }
    }

    pub fn next_int(&mut self) -> i32 {
        self.last_num = self.last_num.wrapping_mul(self.a).wrapping_add(self.c); // modulus is treated as 2^32, 32 bit overflow will automatically wrap around
        self.last_num & 0x7FFFFFFF // ignore bit sign bit to ensure non-negative output
    }

    pub fn next_double(&mut self) -> f64 {
        self.next_int() as f64 / 4294967296.0 // 2^32
    }

   
}
