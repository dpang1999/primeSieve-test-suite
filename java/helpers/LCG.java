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
        lastNumber = (a * lastNumber + c) % m;
        return lastNumber;
    }

    public double nextDouble() {
        return (double) nextInt() / m;
    }
}
