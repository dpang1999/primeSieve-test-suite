// Linear Congruential Generator (LCG)
export class LCG {
    private last_num: number;
    private m: number;
    private a: number;
    private c: number;

    constructor(seed: number, m: number, a: number, c: number) {
        this.last_num = seed;
        this.m = m;
        this.a = a;
        this.c = c;
    }

    nextDouble(): number {
        return this.nextInt() / this.m;
    }

    nextInt(): number {
        this.last_num = (this.a * this.last_num + this.c) % this.m;
        return this.last_num;
    }

}