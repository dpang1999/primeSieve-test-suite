package helpers;
import java.util.ArrayList;

public class LCG {
    private int lastNumber;
    private int m;
    private int a;
    private int c;
    public LCG(int seed, int m, int a, int c) {
        this.m = m;
        this.a = a;
        this.c = c;
        this.lastNumber = seed;
    }

    public int nextInt() {
        lastNumber = (a * lastNumber + c); // modulus is treated as 2^32, 32 bit overflow will automatically wrap around
        return lastNumber & 0x7FFFFFFF; // ignore bit sign bit to ensure non-negative output
    }

    public double nextDouble() {
        return (double) nextInt() / 4294967296.0; // 2^32
    }
}

