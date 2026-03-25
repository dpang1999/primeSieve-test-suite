// Linear Congruential Generator (LCG)
export class LCG {
    constructor(seed, m, a, c) {
        this.last_num = seed;
        this.m = m;
        this.a = a;
        this.c = c;
    }
    nextDouble() {
        return this.nextInt() / this.m;
    }
    nextInt() {
        // modulus should be a very slow operation as TypeScript/JavaScript does not have native support for it
        this.last_num = (this.a * this.last_num + this.c) % this.m;
        return this.last_num;
    }
}
