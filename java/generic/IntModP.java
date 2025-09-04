package generic;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Formatter;
import java.util.List;

public class IntModP implements IField<IntModP>,
		IOrdered<IntModP>, ICopiable<IntModP>, IPrimitiveRoots<IntModP>, IMath<IntModP> {
	long d;
    long p;

	public boolean printShort = true;

	//public static int fCount;
    IntModP(long d, long p) {
            this.d = d % p;
            this.p = p;
    }
    
    public IntModP copy() {
        return new IntModP(d, p);
    }
    public IntModP[] newArray(int size) {
        return new IntModP[size];
    }
    public IntModP a(IntModP o) {
        //fCount++;
        if (o == null)
            return new IntModP(d, p);
        if (Long.MAX_VALUE - d < o.d) {
            //System.out.println("OVERFLOW: d=" + d + ", o.d=" + o.d);
            BigInteger bigD = BigInteger.valueOf(d);
            BigInteger bigOD = BigInteger.valueOf(o.d);
            BigInteger bigP = BigInteger.valueOf(p);
            long result = bigD.add(bigOD).mod(bigP).longValue();
            return new IntModP(result, p);
        }
        else
            return new IntModP((d + o.d) % p, p);
    }
    public void ae(IntModP o) {
        //fCount++;
        if (o != null)
            d = (d + o.d) % p;
        if (Long.MAX_VALUE - d < o.d) {
            //System.out.println("OVERFLOW: d=" + d + ", o.d=" + o.d);
            BigInteger bigD = BigInteger.valueOf(d);
            BigInteger bigOD = BigInteger.valueOf(o.d);
            BigInteger bigP = BigInteger.valueOf(p);
            long result = bigD.add(bigOD).mod(bigP).longValue();
            d = result;
        }
    }
    public IntModP s(IntModP o) {
        //fCount++;
        if (o != null)
            return new IntModP((d - o.d + p) % p, p);   
        else
            return new IntModP(d, p);
    }
    public void se(IntModP o) {
        //fCount++;
        if (o != null)
            d = (d - o.d + p) % p;
    }
    public IntModP m(IntModP o) {
        //fCount++;
        if (o == null) {
            return new IntModP(0, p);
        }
        if (d != 0 && (Long.MAX_VALUE / d) + 1 < o.d) { // Add 1 to compensate for edge case with integer division
            BigInteger bigD = BigInteger.valueOf(d);
            BigInteger bigOD = BigInteger.valueOf(o.d);
            BigInteger bigP = BigInteger.valueOf(p);
            long result = bigD.multiply(bigOD).mod(bigP).longValue();
            return new IntModP(result, p);
        }
        return new IntModP((d * o.d) % p, p);

            
    }
    public void me(IntModP o) {
        //fCount++;
        if (d != 0 && (Long.MAX_VALUE / d) + 1 < o.d) { // Add 1 to compensate for edge case with integer division
            BigInteger bigD = BigInteger.valueOf(d);
            BigInteger bigOD = BigInteger.valueOf(o.d);
            BigInteger bigP = BigInteger.valueOf(p);
            d = bigD.multiply(bigOD).mod(bigP).longValue();
        }
        else 
            d = (d * o.d) % p;
    }

    public IntModP d(IntModP o) {
        if (o == null || o.d == 0) {
            throw new ArithmeticException("Division by zero in IntModP");
        }
        long inverse = modInverse(o.d, p); // Compute modular inverse of o.d
        if (d != 0 && (Long.MAX_VALUE / d) + 1 < inverse) // Add 1 to compensate for edge case with integer division
        {
            System.out.println("OVERFLOW: d=" + d + ", inverse=" + inverse);
            BigInteger bigD = BigInteger.valueOf(d);
            BigInteger bigInverse = BigInteger.valueOf(inverse);
            BigInteger bigP = BigInteger.valueOf(p);
            d = bigD.multiply(bigInverse).mod(bigP).longValue();
            return new IntModP(d, p);
        }
        return new IntModP((d * inverse) % p, p); // Multiply by the inverse modulo p
    }
    
    public void de(IntModP o) {
        if (o == null || o.d == 0) {
            throw new ArithmeticException("Division by zero in IntModP");
        }
        long inverse = modInverse(o.d, p); // Compute modular inverse of o.d
        long temp = d;
        d = (d * inverse) % p; // Update the current value
        if (d != 0 && (Long.MAX_VALUE / temp) + 1 < inverse) // Add 1 to compensate for edge case with integer division
        {
            BigInteger bigD = BigInteger.valueOf(d);
            BigInteger bigInverse = BigInteger.valueOf(modInverse(o.d, p));
            BigInteger bigP = BigInteger.valueOf(p);

            // Compute (d * inverse) % p using BigInteger
            d = bigD.multiply(bigInverse).mod(bigP).longValue();
        }
    }

    private static long modInverse(long b, long p) {
        long t = 0, newT = 1;
        long r = p, newR = b;

        while (newR != 0) {
            long quotient = r / newR;

            // Update t and r
            long tempT = t;
            t = newT;
            newT = tempT - quotient * newT;

            long tempR = r;
            r = newR;
            newR = tempR - quotient * newR;
        }

        if (r > 1) {
            throw new ArithmeticException("b is not invertible modulo p");
        }
        if (t < 0) {
            t += p;
        }

        return t;
    }


    public IntModP coerce(int i) {
        return new IntModP((i % p + p) % p, p);
    }

    public IntModP coerce(double i) {
        return new IntModP(((int) i % p + p) % p, p);
    }
    public double coerce() {
        return d;
    }
    public boolean isZero() {
        return d == 0;
    }
    public boolean isOne() {
        return d == 1;
    }
    public IntModP zero() {
        return new IntModP(0, p);
    }
    public IntModP one() {
        return new IntModP(1, p);
    }

    public IntModP primitiveRoot(long n) {
        if (n <= 0 || n >= p) {
            throw new IllegalArgumentException("n must be in range [1, p-1]");
        }

        // Factorize p-1
        List<Long> factors = factorize(p - 1);

        // Test candidates for primitive root
        for (long g = 2; g < p; g++) {
            boolean isRoot = true;
            for (long factor : factors) {
                if (modPow(g, (p - 1) / factor, p) == 1) {
                    isRoot = false;
                    break;
                }
            }
            if (isRoot) {
                return new IntModP(g, p);
            }
        }
        throw new IllegalArgumentException("No primitive root found");

    }
    private List<Long> factorize(long n) {
        List<Long> factors = new ArrayList<>();
        for (long i = 2; i * i <= n; i++) {
            if (n % i == 0) {
                factors.add(i);
                while (n % i == 0) {
                    n /= i;
                }
            }
        }
        if (n > 1) {
            factors.add(n);
        }
        return factors;
    }

    public IntModP pow(long exp) {
        if (exp < 0) {
            throw new IllegalArgumentException("Exponent must be non-negative");
        }
        return new IntModP(modPow(d, exp, p), p);
    }

    public IntModP[] precomputeRootsOfUnity(int n, int direction) {
		// Ensure n divides (p - 1)
		if ((p - 1) % n != 0) {
            throw new IllegalArgumentException("n must divide p-1 for roots of unity to exist");
		}

		// Find a primitive root modulo p
		IntModP primitiveRoot = primitiveRoot(p - 1);
        //System.out.println("Primitive root: " + primitiveRoot);

		// Compute the primitive n-th root of unity
		IntModP omega = primitiveRoot.pow((p - 1) / n);
        //System.out.println("Omega: " + omega + " p: " + p + " n: " + n + " p-1/n= " + (p - 1) / n);

		// Generate all n-th roots of unity
		IntModP[] roots = new IntModP[n];
		for (int k = 0; k < n; k++) {
			// Compute omega^k * direction
			long exponent = (k * direction) % (p - 1);
			if (exponent < 0) {
				exponent += (p - 1); // Ensure positive exponent
			}
			roots[k] = omega.pow(exponent);
            //System.out.println("Root " + k + ": " + roots[k] + " exponent: " + exponent);
		}

		return roots;
	}

    public int compareTo(IntModP o) {
        if (o == null)
            return 1;
        if (d < o.d)
            return -1;
        else if (d > o.d)
            return 1;
        else
            return 0;
    }
    public String toString() {
        if (printShort) {
            try (Formatter fmt = new Formatter()) {
                fmt.format("%6d", d);
                return fmt.toString();
            }
        } else {
            return Long.toString(d);
        }
    }
    public boolean lt(IntModP o) {
        if (o == null)
            return false;
        return d < o.d;
    }
    public boolean le(IntModP o) {
        if (o == null)
            return false;
        return d <= o.d;  
    }
    public boolean gt(IntModP o) {
        if (o == null)
            return false;
        return d > o.d;
    }
    public boolean ge(IntModP o) {
        if (o == null)
            return false;
        return d >= o.d;
    }

    public boolean eq(IntModP o) {
        if (o == null)
            return false;
        return d == o.d && p == o.p;
    }

    private static long modPow(long base, long exp, long mod) {
        if (mod <= 0) {
            throw new IllegalArgumentException("Modulus must be positive");
        }
        long result = 1;
        base = base % mod; // Ensure base is within the range of mod
        if (Long.MAX_VALUE / base  < base) {
            //System.out.println("Too big, do BigInteger modPow instead");
            BigInteger bigBase = BigInteger.valueOf(base);
            BigInteger bigMod = BigInteger.valueOf(mod);
            BigInteger bigExp = BigInteger.valueOf(exp);

            // Use BigInteger's modPow method
            return bigBase.modPow(bigExp, bigMod).longValue();
        }
        while (exp > 0) {
            // If exp is odd, multiply result with base
            if ((exp & 1) == 1) {
                long temp = result;
                result = (result * base) % mod;
                // Check for overflow
                if (base != 0 && Long.MAX_VALUE / base < temp) {
                    //System.out.println("OVERFLOW: result=" + temp + ", base=" + base + ", mod=" + mod);
                    BigInteger bigTemp = BigInteger.valueOf(temp);
                    BigInteger bigBase = BigInteger.valueOf(base);
                    BigInteger bigResult = bigTemp.multiply(bigBase).mod(BigInteger.valueOf(mod));
                    result = bigResult.longValue();
                }
            }

            // Square the base and reduce modulo mod
            long temp = base;
            base = (base * base) % mod;
            // Check for overflow
            if (temp != 0 && Long.MAX_VALUE / temp < temp) {
                //System.out.println("OVERFLOW: base=" + temp + ", exp=" + exp + ", mod=" + mod);
                BigInteger bigTemp = BigInteger.valueOf(temp);
                BigInteger bigMod = BigInteger.valueOf(mod);
                BigInteger bigResult = bigTemp.multiply(bigTemp).mod(bigMod);
                base = bigResult.longValue();
            }

            // Divide exp by 2
            exp >>= 1;
        }

        return result;
    }

    public void sqrt() {
        if (d == 0) return; // 0 is a square root of 0
        // Use Tonelli-Shanks algorithm or similar to find square root
        throw new UnsupportedOperationException("Square root not implemented");
    }

    public void abs() {
        if (d < 0) {
            d = -d;
        }
    }
}