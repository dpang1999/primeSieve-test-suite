package generic;
import java.util.*;
import java.lang.reflect.*;

/**
 * Computes FFT's of complex, double precision data where n is an integer power
 * of 2. This appears to be slower than the Radix2 method, but the code is
 * smaller and simpler, and it requires no extra storage.
 * <P>
 * 
 * @author Laurentiu Dragan ldragan@scl.csd.uwo.ca,
 * @author Bruce R. Miller bruce.miller@nist.gov,
 * @author Derived from GSL (Gnu Scientific Library),
 * @author GSL's FFT Code by Brian Gough bjg@vvv.lanl.gov
 */

/*
 * See {@link ComplexDoubleFFT ComplexDoubleFFT} for details of data layout.
 */

public class GenFFT<N extends IField<N> & IOrdered<N> & ICopiable<N> & IPrimitiveRoots<N>> {

	private N c;

	final static boolean DEBUG = false;
	final static boolean MUTABLE = false;

	/*GenFFT(N num, ComplexField<N> c) {
		this.num = num;
		this.c = c;
	}*/
	GenFFT (N data)
	{
		this.c = data;
	}

	public static final double num_flops(int N) {
		double Nd = (double) N;
		double logN = (double) log2(N);

		return (5.0 * Nd - 2) * logN + 2 * (Nd + 1);
	}

	/** Compute Fast Fourier Transform of (complex) data, in place. */
	public void transform(N[] data) {
		transform_internal(data, -1);
	}

	/* Compute Inverse FFT of (complex) data, in place. */
	public void inverse(N[] data) {
		transform_internal(data, +1);
		// Normalize
		int nd = data.length;
		N test = c.coerce(nd);
		for (int i = 0; i < nd; i++)
			data[i].de(test);
	}

	/**
	 * Accuracy check on FFT of data. Make a copy of data, Compute the FFT, then
	 * the inverse and compare to the original. Returns the rms difference.
	 */
	public double test(N[] data) {
		int nd = data.length;
		// Make duplicate for comparison
		ArrayList<N> copy = new ArrayList<N>(nd);
		for (int i = 0; i < nd; i++) {
			copy.add(data[i].copy());
		}
		// double copy[] = new double[nd];
		// System.arraycopy(data,0,copy,0,nd);
		// Transform & invert
		transform(data);

		System.out.println("After transform:" + Arrays.toString(data));

		inverse(data);

		System.out.println("After inverse:" + Arrays.toString(data));
		// Compute RMS difference. TODO or not idk
		double diff = 0.0;
		/*for (int i = 0; i < nd; i++) {
			N d = data[i].copy();
			double real = d.re.coerce();
			double imag = d.im.coerce();
			double realDiff = real - copy.get(i).re.coerce();
			double imagDiff = imag - copy.get(i).im.coerce();
			diff += realDiff * realDiff + imagDiff * imagDiff;
		}*/

		return Math.sqrt(diff/(nd*2)); // nd*2 as the original had double array size for imaginary and complex components separately

		
	}

	/** Make a random array of n (complex) elements. */
	@SuppressWarnings("unchecked")
	public ComplexField<DoubleField>[] makeRandom(int n) {
		Random random = new Random(12345);
		ComplexField<DoubleField>[] data = (ComplexField<DoubleField>[]) Array.newInstance(c.getClass(), n);
		DoubleField temp = new DoubleField(0);
		for (int i = 0; i < n; i++) {
			data[i] = new ComplexField(temp.coerce(random.nextDouble()),
					temp.coerce(random.nextDouble()));
		}
		return data;
	}

	/** Simple Test routine. */
	public static void main(String args[]) {
		ComplexField<DoubleField> c = new ComplexField<DoubleField>(new DoubleField(0),
				new DoubleField(0));
		DoubleField num = new DoubleField(0);
		GenFFT<ComplexField<DoubleField>> fft = new GenFFT<ComplexField<DoubleField>>(
				c);
		int type = 0;
		if (args.length == 0 && type == 0) {
			int n = 16;
			ComplexField<DoubleField>[] data = fft.makeRandom(n);

			/*ComplexField<DoubleField>[] data2 = (ComplexField<DoubleField>[]) java.lang.reflect.Array.newInstance(ComplexField.class, 4);
			data2[0] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));
			data2[1] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));
			data2[2] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));
			data2[3] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));*/
			System.out.println(Arrays.toString(data));
			System.out.println("n=" + n + " => RMS Error="+ fft.test(data));


			double[] in1 = {38.0,  0.0, 44.0, 87.0,  6.0, 45.0, 22.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			ComplexField<DoubleField>[] data1 = (ComplexField<DoubleField>[]) Array.newInstance(c.getClass(), in1.length / 2);
      		double[] in2 = {80.0, 18.0, 62.0, 90.0, 17.0, 96.0, 27.0, 97.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			ComplexField<DoubleField>[] data2 = (ComplexField<DoubleField>[]) Array.newInstance(c.getClass(), in2.length / 2);
			for(int i = 0; i < in1.length; i += 2) {
				data1[i / 2] = new ComplexField<DoubleField>(
					new DoubleField(in1[i]), new DoubleField(in1[i + 1]));
				data2[i / 2] = new ComplexField<DoubleField>(
					new DoubleField(in2[i]), new DoubleField(in2[i + 1]));
			}
			System.out.println("Data1: " + Arrays.toString(data1));
			System.out.println("Data2: " + Arrays.toString(data2));
			fft.transform(data1);
			System.out.println("Transformed Data1: " + Arrays.toString(data1));
			fft.transform(data2);
			System.out.println("Transformed Data2: " + Arrays.toString(data2));

			ComplexField<DoubleField>[] product = (ComplexField<DoubleField>[]) Array.newInstance(c.getClass(), data1.length);
			// multiply the complex numbers
			for (int i = 0; i < data1.length; i++) {
				product[i] = data1[i].m(data2[i]);
			}
			System.out.println("Product: " + Arrays.toString(product));
			fft.inverse(product);
			System.out.println("Inverse Product: " + Arrays.toString(product));
			fft.inverse(data1);
			System.out.println("Inverse Data1: " + Arrays.toString(data1));
			fft.inverse(data2);
			System.out.println("Inverse Data2: " + Arrays.toString(data2));
		}
		else if (type == 1) {
			int[] basic = {1,2,3,4};
			IntModP basicNum = new IntModP(1, 13);
			GenFFT basicFft = new GenFFT(basicNum);
			IntModP[] basicData = new IntModP[basic.length];
			for(int i = 0; i < basic.length; i ++) {
				basicData[i] = basicNum.coerce(basic[i]);
			}
			System.out.println("Basic Data: " + Arrays.toString(basicData));
			basicFft.transform(basicData);
			System.out.println("Transformed Basic Data: " + Arrays.toString(basicData));
			basicFft.inverse(basicData);
			System.out.println("Inverse Basic Data: " + Arrays.toString(basicData));


			int[] in1 = {38,  0, 44, 87,  6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0};
			int[] in2 = {80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0};
			int[] out = { 3040,    684,  5876, 11172,  5420, 16710, 12546, 20555,
                      16730, 15704, 21665,  5490, 13887,  4645,  9021,0 };
			IntModP finiteNum = new IntModP(1, 40961);
			GenFFT finiteFft = new GenFFT(finiteNum);
			IntModP[] data1 = new IntModP[in1.length];
			IntModP[] data2 = new IntModP[in2.length];
			for(int i = 0; i < in1.length; i ++) {
				data1[i] = finiteNum.coerce(in1[i]);
				data2[i] = finiteNum.coerce(in2[i]);
			}
			System.out.println("Data1: " + Arrays.toString(data1));
			System.out.println("Data2: " + Arrays.toString(data2));
			finiteFft.transform(data1);
			System.out.println("Transformed Data1: " + Arrays.toString(data1));
			finiteFft.transform(data2);
			System.out.println("Transformed Data2: " + Arrays.toString(data2));

			IntModP[] product = new IntModP[data1.length];
			// multiply the finite field numbers
			for (int i = 0; i < data1.length; i++) {
				product[i] = data1[i].m(data2[i]);
			}
			System.out.println("Product: " + Arrays.toString(product));
			finiteFft.inverse(product);
			System.out.println("Inverse Product: " + Arrays.toString(product));
			finiteFft.inverse(data1);
			System.out.println("Inverse Data1: " + Arrays.toString(data1));
			finiteFft.inverse(data2);
			System.out.println("Inverse Data2: " + Arrays.toString(data2));
		}
		for (int i = 0; i < args.length; i++) {
			int n = Integer.parseInt(args[i]);
			ComplexField<DoubleField>[] data = fft.makeRandom(n);
			System.out.println(Arrays.toString(data));
			System.out.println("n=" + n + " => RMS Error="
					+ fft.test(data));
		}
	}

	/* ______________________________________________________________________ */

	protected static int log2(int n) {
		int log = 0;
		for (int k = 1; k < n; k *= 2, log++)
			;
		if (n != (1 << log))
			throw new Error("FFT: Data length is not a power of 2!: " + n);
		return log;
	}

	
	/*private ComplexField<N>[] precomputeRootsOfUnity(int n, int direction) {
		ComplexField<N>[] roots = (ComplexField<N>[]) Array.newInstance(c.getClass(), n);
		for (int k = 0; k < n; k++) {
			double angle = 2.0 * Math.PI * k * direction / n;
			N realPart = c.re.coerce(Math.cos(angle));
			N imagPart = c.re.coerce(Math.sin(angle));
			roots[k] = new ComplexField<>(realPart, imagPart);
		}
		return roots;
	}
	private IntModP[] precomputeRootsOfUnity(int n, int direction, IntModP num) {
		// Ensure n divides (p - 1)
		if ((num.p - 1) % n != 0) {
			throw new IllegalArgumentException("n must divide p-1 for roots of unity to exist");
		}

		// Find a primitive root modulo p
		IntModP primitiveRoot = num.primitiveRoot(num.p - 1);

		// Compute the primitive n-th root of unity
		IntModP omega = primitiveRoot.pow((num.p - 1) / n);

		// Generate all n-th roots of unity
		IntModP[] roots = new IntModP[n];
		for (int k = 0; k < n; k++) {
			// Compute omega^k * direction
			int exponent = (k * direction) % (num.p - 1);
			if (exponent < 0) {
				exponent += (num.p - 1); // Ensure positive exponent
			}
			roots[k] = omega.pow(exponent);
		}

		return roots;
	}*/


	protected void transform_internal(N[] data, int direction) {
		if (data.length == 0)
			return;
		int n = data.length;
		if (n == 1)
			return; // Identity operation
		int logn = log2(n);

		// bit reverse
		bitreverse(data);

		// apply fft recursion
		// this loop is executed log2(N) times
		N w;
		int dual = 1;
		N[] roots = data[0].precomputeRootsOfUnity(n, direction);
		//System.out.println("Roots: " + Arrays.toString(roots));
		for (int bit = 0; bit < logn; bit++, dual *= 2) {
			for (int a = 0; a < dual; a++) {
				w = roots[a * (n / (2 * dual))]; // Use precomputed root
				for (int b = 0; b < n; b += 2 * dual) {
					int i = b + a;
					int j = b + a + dual;

					// Twiddle factor multiplication
					N z1 = data[j].m(w);

					// Butterfly operation
					N tempI = data[i];
					data[i] = tempI.a(z1);   // Add
					data[j] = tempI.s(z1);   // Subtract
				}
				//System.out.println("After a=" + a + ": " + Arrays.toString(data));
			}
		}


	}
/* 
	protected void transform_internal(ComplexField<N>[] data, int direction) {
		if (data.length == 0)
			return;
		int n = data.length;
		if (n == 1)
			return; // Identity operation!
		int logn = log2(n);

		// bit reverse the input data for decimation in time algorithm 
		bitreverse(data);

		// apply fft recursion 
		// this loop executed log2(N) times 
		final N n1 = c.re.coerce(1);
		final N n0 = c.re.coerce(0);
		final ComplexField<N> c10 = new ComplexField<N>(c.re.coerce(1), c.re.coerce(0));
		final N n2 = c.re.coerce(2);
		N theta, s, t, s2;
		ComplexField<N> w;
		int dual = 1;
		ComplexField<N>[] roots = precomputeRootsOfUnity(n, direction);
		for (int bit = 0; bit < logn; bit++, dual *= 2) {
			// Roots of unity variant
			for (int a = 0; a < dual; a++) {
				w = roots[a * (n / (2 * dual))]; // Use precomputed root
				for (int b = 0; b < n; b += 2 * dual) {
					int i = b + a;
					int j = b + a + dual;

					ComplexField<N> z1 = data[j].m(w); // Twiddle factor multiplication
					data[j] = data[i].s(z1);          // Subtract
					data[i] = data[i].a(z1);          // Add
				}
			}
		}
		

	}
	protected void transform_internal(IntModP[] data, int direction) {
		if (data.length == 0)
			return;
		int n = data.length;
		if (n == 1)
			return; // Identity operation!
		int logn = log2(n);

		// bit reverse the input data for decimation in time algorithm
		bitreverse(data);
		//System.out.println("After reverse:" + Arrays.toString(data));

		// apply fft recursion 
		IntModP w;
		int dual = 1;
		IntModP[] roots = precomputeRootsOfUnity(n, direction, num);
		//System.out.println("Roots: " + Arrays.toString(roots));
		for (int bit = 0; bit < logn; bit++, dual *= 2) {
			for (int a = 0; a < dual; a++) {
				w = roots[a * (n / (2 * dual))]; // Use precomputed root
				for (int b = 0; b < n; b += 2 * dual) {
					int i = b + a;
					int j = b + a + dual;

					// Twiddle factor multiplication
					IntModP z1 = data[j].m(w);

					// Butterfly operation
					IntModP tempI = data[i];
					data[i] = tempI.a(z1);   // Add
					data[j] = tempI.s(z1);   // Subtract
				}
				//System.out.println("After a=" + a + ": " + Arrays.toString(data));
			}
		}
			
	}
	 */
	/*
	protected void bitreverse(ComplexField<N>[] data) {
		// This is the Goldrader bit-reversal algorithm 
		int n = data.length;
		int nm1 = n - 1;
		int i = 0;
		int j = 0;
		for (; i < nm1; i++) {

			int k = n >> 1;

			if (i < j) {
				ComplexField<N> tmp = data[i];
				data[i] = data[j];
				data[j] = tmp;
			}

			while (k <= j) {
				// j = j - k ;
				j -= k;

				// k = k / 2 ;
				k >>= 1;
			}
			j += k;
		}
	}*/

	protected void bitreverse(N[] data) {
		/* This is the Goldrader bit-reversal algorithm */
		int n = data.length;
		int nm1 = n - 1;
		int i = 0;
		int j = 0;
		for (; i < nm1; i++) {

			if (i < j) {
				N tmp = data[i];
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
