package generic;

public class BitPackedExponent implements IExponent<BitPackedExponent> {
    // Stores degree in top 16 bits, exponents in lower 48 bits (6 x 8 bits)
    public long packed;

    public BitPackedExponent(long packed) {
        this.packed = packed;
    }

    // Construct from array of 6 exponents
    public static BitPackedExponent fromArray(int[] exponents) {
        long packed = 0;
        int degree = 0;
        for (int i = 0; i < 6; i++) {
            int exp = exponents[i] & 0xFF;
            int shift = 40 - (8 * i);
            packed |= ((long) exp) << shift;
            degree += exp;
        }
        packed |= ((long) degree & 0xFFFF) << 48;
        return new BitPackedExponent(packed);
    }

    // Get exponents as array
    public int[] getExponents() {
        int[] exps = new int[6];
        for (int i = 0; i < 6; i++) {
            int shift = 40 - (8 * i);
            exps[i] = (int) ((packed >> shift) & 0xFF);
        }
        return exps;
    }

    public BitPackedExponent add(BitPackedExponent o) {
        return new BitPackedExponent(this.packed + o.packed);
    }

    public BitPackedExponent sub(BitPackedExponent o) {
        return new BitPackedExponent(this.packed - o.packed);
    }

    public int degree() {
        return (int) ((packed >> 48) & 0xFFFF);
    }

    // LCM
    public BitPackedExponent lcm(BitPackedExponent o) {
        long selfExponents = this.packed & 0x0000FFFFFFFFFFFFL;
        long otherExponents = o.packed & 0x0000FFFFFFFFFFFFL;
        long lcmExponents = 0;
        int degree = 0;
        for (int i = 0; i < 6; i++) {
            int shift = 40 - (8 * i);
            int a = (int) ((selfExponents >> shift) & 0xFF);
            int b = (int) ((otherExponents >> shift) & 0xFF);
            int maxExp = Math.max(a, b);
            lcmExponents |= ((long) maxExp) << shift;
            degree += maxExp;
        }
        lcmExponents |= ((long) degree & 0xFFFF) << 48;
        return new BitPackedExponent(lcmExponents);
    }

    // Lexicographic compare
    public int compareTo(BitPackedExponent o) {
        int[] a = this.getExponents();
        int[] b = o.getExponents();
        for (int i = 0; i < 6; i++) {
            if (a[i] < b[i]) return -1;
            if (a[i] > b[i]) return 1;
        }
        return 0;
    }

    // Can reduce
    public boolean canReduce(BitPackedExponent o) {
        int[] a = this.getExponents();
        int[] b = o.getExponents();
        for (int i = 0; i < 6; i++) {
            if (a[i] < b[i]) return false;
        }
        return true;
    }

    // Display
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Degree: ").append(degree()).append(", Exponents: ");
        int[] exps = getExponents();
        for (int e : exps) sb.append(String.format("%02X ", e));
        return sb.toString();
    }
}