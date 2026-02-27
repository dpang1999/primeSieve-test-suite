// Linear Congruential Generator (LCG)
export class LCG {
    private last_num: number;
    private m: number;
    private a: number;
    private c: number;

    constructor(seed: number, m: number, a: number, c: number) {
        this.last_num = seed;
        this.m = m; 
        // For the sake of fairness, modulus will always be 2^32 to bypass modulus bias between languages
        this.a = a;
        this.c = c;
    }

    nextDouble(): number {
        return this.nextUint32() / 4294967296.0; // 2^32
    }

    nextUint32(): number {
        this.last_num = (this.a * this.last_num + this.c) | 0 // modulus is treated as 2^32, 32 bit overflow will automatically wrap around
        return this.last_num & 0x7FFFFFFF; // ignore bit sign bit to ensure non-negative output
    }

}