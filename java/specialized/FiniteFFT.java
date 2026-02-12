package specialized;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import generic.IntModP;
import helpers.LCG;
import helpers.FindPrime;

/** Computes FFT's of complex, double precision data where n is an integer power of 2.
  * This appears to be slower than the Radix2 method,
  * but the code is smaller and simpler, and it requires no extra storage.
  * <P>
  *
  * @author Bruce R. Miller bruce.miller@nist.gov,
  * @author Derived from GSL (Gnu Scientific Library), 
  * @author GSL's FFT Code by Brian Gough bjg@vvv.lanl.gov
  */

  /* See {@link ComplexDoubleFFT ComplexDoubleFFT} for details of data layout.
   */

public class FiniteFFT {

  public static final double num_flops(int N)
  {
	 double Nd = (double) N;
	 double logN = (double) log2(N);

	 return (5.0*Nd-2)*logN + 2*(Nd+1);
   }
    private static int modInverse(int b, int p) {
        int t = 0, newT = 1;
        int r = p, newR = b;

        while (newR != 0) {
            int quotient = r / newR;

            // Update t and r
            int tempT = t;
            t = newT;
            newT = tempT - quotient * newT;

            int tempR = r;
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

    public static int primitiveRoot(int modulus) {
        List<Integer> factors = factorize(modulus - 1);

        for (int g = 2; g < modulus; g++) {
            boolean isRoot = true;
            for (int factor : factors) {
                if (modPow(g, (modulus - 1) / factor, modulus) == 1) {
                    isRoot = false;
                    break;
                }
            }
            if (isRoot) {
                return g;
            }
        }
        throw new IllegalArgumentException("No primitive root found");

    }

    private static List<Integer> factorize(int n) {
        List<Integer> factors = new ArrayList<>();
        for (int i = 2; i * i <= n; i++) {
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


    private static int modPow(int base, int exp, int mod) {
        if (mod <= 0) {
            throw new IllegalArgumentException("Modulus must be positive");
        }
        int result = 1;
        base = base % mod; // Ensure base is within the range of mod
        while (exp > 0) {
            // If exp is odd, multiply result with base
            if ((exp & 1) == 1) {
                long temp = result;
                result = (result * base) % mod;
            }

            // Square the base and reduce modulo mod
            base = (base * base) % mod;
            // Divide exp by 2
            exp >>= 1;
        }

        return result;
    }

    public static int[] precomputeRootsOfUnity(int n, int direction, int modulus) {
		// Ensure n divides (p - 1)
		if ((modulus - 1) % n != 0) {
            throw new IllegalArgumentException("n must divide p-1 for roots of unity to exist");
		}

		// Find a primitive root modulo p
		int primitiveRoot = primitiveRoot(modulus);
        //System.out.println("Primitive root: " + primitiveRoot);

		// Compute the primitive n-th root of unity
		int omega = modPow(primitiveRoot, (modulus - 1) / n, modulus);
        //System.out.println("Omega: " + omega + " p: " + p + " n: " + n + " p-1/n= " + (p - 1) / n);

		// Generate all n-th roots of unity
		int[] roots = new int[n];
		for (int k = 0; k < n; k++) {
			// Compute omega^k * direction
			int exponent = (k * direction) % (modulus - 1);
			if (exponent < 0) {
				exponent += (modulus - 1); // Ensure positive exponent
			}
			roots[k] = modPow(omega, exponent, modulus);
            //System.out.println("Root " + k + ": " + roots[k] + " exponent: " + exponent);
		}

		return roots;
	}

  /** Compute Fast Fourier Transform of (complex) data, in place.*/
  public static void transform (long data[], int modulus, int root) {
    transform_internal(data, -1, modulus, root); }

  /** Compute Inverse Fast Fourier Transform of (complex) data, in place.*/
  public static void inverse (long data[], int modulus, int root) {
    transform_internal(data, +1, modulus, root);
    // Normalize
    int nd=data.length;
    int n =nd;
    for(int i=0; i<nd; i++)
      data[i] = data[i] * modInverse(n, modulus) % modulus;
  }

   /** Accuracy check on FFT of data. Make a copy of data, Compute the FFT, then
    * the inverse and compare to the original.  Returns the rms difference.*/
  /*public static double test(long data[]){
    int nd = data.length;
    // Make duplicate for comparison
    double copy[] = new double[nd];
    System.arraycopy(data,0,copy,0,nd);
    // Transform & invert
    transform(data);
    System.out.println("After transform:" + Arrays.toString(data));
    inverse(data);
    System.out.println("After inverse:" + Arrays.toString(data));
    // Compute RMS difference.
    double diff = 0.0;
    for(int i=0; i<nd; i++) {
      double d = data[i]-copy[i];
      diff += d*d; 
    }

    return Math.sqrt(diff/nd); } */

  /** Make a random array of n (complex) elements. */
  public static double[] makeRandom(int n){
    LCG random = new LCG(12345, 1345, 16645, 1013904);
    int nd = 2*n;
    double data[] = new double[nd];
    for(int i=0; i<nd; i++)
      data[i]= random.nextDouble();
    return data; }

  /** Simple Test routine. */
  public static void main(String args[]){
    int mode = 1;
    if (mode == 0) {
      int n = Integer.parseInt(args[0]);
      System.out.println("FFT with " + n);
      int modulus = FindPrime.findPrimeCongruentOneModN(n);
      int root = primitiveRoot(modulus);
      LCG rand = new LCG(12345, 1345, 16645, 1013904);
      long data[] = new long[n];
      for (int i = 0; i < n; i=i+2) {
          long real = (long) rand.nextInt();
          data[i] = real;
          }
      for (int loop = 1; loop <= 10; loop++) { //looping for JIT warmup

        transform(data, modulus, root);
        inverse(data, modulus, root);
        System.out.println("Loop " + loop + " done");
      }
    }
    else {

        long[] in1 = {38L,  0L, 44L, 87L,  6L, 45L, 22L, 93L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L};
        long[] in2 = {80L, 18L, 62L, 90L, 17L, 96L, 27L, 97L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L};
        int prime = 40961;
        int root = primitiveRoot(prime);
        transform(in1, prime, root);
        System.out.println("Transformed in1: " + Arrays.toString(in1));
        transform(in2, prime, root);
        System.out.println("Transformed in2: " + Arrays.toString(in2));

        long[] product = new long[in1.length];
        // multiply the complex numbers
        for (int i = 0; i < in1.length; i ++) {
            product[i] = in1[i] * in2[i] % prime;
        }
        System.out.println("Product: " + Arrays.toString(product));
        inverse(product, prime, root);
        System.out.println("Inverse Product: " + Arrays.toString(product));
        inverse(in1, prime, root);
        System.out.println("Inverse in1: " + Arrays.toString(in1));
        inverse(in2, prime, root);
        System.out.println("Inverse in2: " + Arrays.toString(in2));
        }
    }


  /* ______________________________________________________________________ */

  protected static int log2 (int n){
    int log = 0;
    for(int k=1; k < n; k *= 2, log++);
    if (n != (1 << log))
      throw new Error("FFT: Data length is not a power of 2!: "+n);
    return log; }

  protected static void transform_internal (long data[], int direction, int modulus, int root) {
	if (data.length == 0) return;    
	int n = data.length;
    if (n == 1) return;         // Identity operation!
    int logn = log2(n);

    /* bit reverse the input data for decimation in time algorithm */
    bitreverse(data) ;
    
    int[] roots = precomputeRootsOfUnity(n, direction, modulus);

    int dual = 1;
    int w;
    for (int bit = 0; bit < logn; bit++, dual *= 2) {
			for (int a = 0; a < dual; a++) {
				w = roots[a * (n / (2 * dual))]; // Use precomputed root
				for (int b = 0; b < n; b += 2 * dual) {
					int i = b + a;
					int j = b + a + dual;

					// Twiddle factor multiplication
					long z1 = data[j]*w % modulus;

					// Butterfly operation
					long tempI = data[i];
					data[i] = (tempI + z1) % modulus;   // Add
					data[j] = (tempI + modulus - z1) % modulus;   // Subtract
				}
			}
    }
  }


  protected static void bitreverse(long data[]) {
    /* This is the Goldrader bit-reversal algorithm */
    int n = data.length;
		int nm1 = n - 1;
		int i = 0;
		int j = 0;
		for (; i < nm1; i++) {

			if (i < j) {
				long tmp = data[i];
				data[i] = data[j];
				data[j] = tmp;
			}

			int k = n >> 1;

			while (k <= j) {
				j -= k;
				k >>= 1;
			}
			j += k;
		}
  }
}








