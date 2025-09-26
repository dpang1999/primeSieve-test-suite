pub struct lcg {
    m: i32,
    a: i32,
    c: i32,
    last_num: i32,
}

impl lcg {
    pub fn new(seed: i32, m: i32, a: i32, c: i32) -> Self {
        Self {
            m,
            a,
            c,
            last_num: seed,
        }
    }

    pub fn next_int(&mut self) -> i32 {
        self.last_num = (self.a * self.last_num + self.c) % self.m;
        self.last_num
    }

    pub fn next_double(&mut self) -> f64 {
        self.next_int() as f64 / self.m as f64
    }

   
}
