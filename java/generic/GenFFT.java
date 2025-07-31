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

public class GenFFT<N extends IField<N> & ITrigonometric<N> & IMath<N> & IOrdered<N> & ICopiable<N>> {

	private ComplexField<N> c;
	//private N num;

	final static boolean DEBUG = false;
	final static boolean MUTABLE = false;

	/*GenFFT(N num, ComplexField<N> c) {
		this.num = num;
		this.c = c;
	}*/
	GenFFT (N re, N im) {
		this.c = new ComplexField<N>(re, im);
	}

	public static final double num_flops(int N) {
		double Nd = (double) N;
		double logN = (double) log2(N);

		return (5.0 * Nd - 2) * logN + 2 * (Nd + 1);
	}

	/** Compute Fast Fourier Transform of (complex) data, in place. */
	public void transform(ComplexField<N>[] data) {
		transform_internal(data, -1);
	}

	/* Compute Inverse FFT of (complex) data, in place. */
	public void inverse(ComplexField<N>[] data) {
		transform_internal(data, +1);
		// Normalize
		int nd = data.length;
		// N norm = num.coerce(1 / (double) nd);
		ComplexField<N> norm = c.coerce(1 / (double) nd);
		for (int i = 0; i < nd; i++)
			data[i].me(norm);
	}

	/**
	 * Accuracy check on FFT of data. Make a copy of data, Compute the FFT, then
	 * the inverse and compare to the original. Returns the rms difference.
	 */
	public double test(ComplexField<N>[] data) {
		int nd = data.length;
		// Make duplicate for comparison
		ArrayList<ComplexField<N>> copy = new ArrayList<ComplexField<N>>(nd);
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
		// Compute RMS difference.
		/*
		 * cannot have static methods in interfaces.
		 */

		/* 
		N diff = num.coerce(0);
		for (int i = 0; i < nd; i++) {
			// C d = data.get(i).substract(copy.get(i));
			C d = data[i];
			d.se(copy.get(i));
			d.re().me(d.re());
			d.im().me(d.im());
			diff.ae(d.re());
			diff.ae(d.im());
		}
		diff.de(num.coerce(nd));
		diff.sqrt();
		return diff;*/


		double diff = 0.0;
		for (int i = 0; i < nd; i++) {
			ComplexField<N> d = data[i].copy();
			double real = d.re.coerce();
			double imag = d.im.coerce();
			double realDiff = real - copy.get(i).re.coerce();
			double imagDiff = imag - copy.get(i).im.coerce();
			diff += realDiff * realDiff + imagDiff * imagDiff;
		}

		return Math.sqrt(diff/(nd*2)); // nd*2 as the original had double array size for imaginary and complex components separately

		
	}

	/** Make a random array of n (complex) elements. */
	@SuppressWarnings("unchecked")
	public ComplexField<N>[] makeRandom(int n) {
		Random random = new Random(12345);
		// int nd = 2*n;
		// ArrayList<C> data = new ArrayList<C>(nd);
		// for(int i=0; i<nd; i++)
		ComplexField<N>[] data = (ComplexField<N>[]) Array.newInstance(c.getClass(), n);
		for (int i = 0; i < n; i++)
			data[i] = new ComplexField<N>(c.re.coerce(random.nextDouble()),
					c.re.coerce(random.nextDouble()));
		return data;
	}

	/** Simple Test routine. */
	public static void main(String args[]) {
		ComplexField<DoubleField> c = new ComplexField<DoubleField>(new DoubleField(0),
				new DoubleField(0));
		DoubleField num = new DoubleField(0);
		GenFFT<DoubleField> fft = new GenFFT<DoubleField>(
				num, num);
		if (args.length == 0) {
			int n = 4;
			ComplexField<DoubleField>[] data = fft.makeRandom(n);

			/*ComplexField<DoubleField>[] data2 = (ComplexField<DoubleField>[]) java.lang.reflect.Array.newInstance(ComplexField.class, 4);
			data2[0] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));
			data2[1] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));
			data2[2] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));
			data2[3] = new ComplexField<DoubleField>(new DoubleField(2), new DoubleField(0));*/
			System.out.println(Arrays.toString(data));
			System.out.println("n=" + n + " => RMS Error="+ fft.test(data));
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

	private ComplexField<N>[] precomputeRootsOfUnity(int n, int direction) {
		ComplexField<N>[] roots = (ComplexField<N>[]) Array.newInstance(c.getClass(), n);
		for (int k = 0; k < n; k++) {
			double angle = 2.0 * Math.PI * k * direction / n;
			N realPart = c.re.coerce(Math.cos(angle));
			N imagPart = c.re.coerce(Math.sin(angle));
			roots[k] = new ComplexField<>(realPart, imagPart);
		}
		return roots;
	}

	protected void transform_internal(ComplexField<N>[] data, int direction) {
		if (data.length == 0)
			return;
		int n = data.length;
		if (n == 1)
			return; // Identity operation!
		int logn = log2(n);

		/* bit reverse the input data for decimation in time algorithm */
		bitreverse(data);

		/* apply fft recursion */
		/* this loop executed log2(N) times */
		final N n1 = c.re.coerce(1);
		final N n0 = c.re.coerce(0);
		final ComplexField<N> c10 = new ComplexField<N>(c.re.coerce(1), c.re.coerce(0));
		final N n2 = c.re.coerce(2);
		N theta, s, t, s2;
		ComplexField<N> w;
		int dual = 1;
		ComplexField<N>[] roots = precomputeRootsOfUnity(n, direction);
		for (int bit = 0; bit < logn; bit++, dual *= 2) {
			/*if (MUTABLE) {
				w = c10.copy();
				theta = c.re.coerce(2.0 * direction * 
					Math.PI / (2.0 * (double) dual));
				s = theta.copy();
				s.sin();
				t = theta.copy();
				t.de(n2);
				t.sin();
				s2 = t.copy();
				s2.me(s2);
				s2.me(n2);
			} else {
				w = c10.copy();
				theta = c.re.coerce(2.0 * direction * 
					Math.PI / (2.0 * (double) dual));
				s = theta.sin();
				//t = theta.d(n2).sin();
				//s2 = t.m(t).m(n2);
				theta.de(n2);
				s2 = theta.sin();
				s2.me(s2);
				s2.me(n2);
			}

			// a = 0 
			for (int b = 0; b < n; b += 2 * dual) {
				int i = b;
				int j = (b + (int) dual);

				if (MUTABLE) {
					ComplexField<N> wd = data[j].copy();

					ComplexField<N> tmp = data[i].copy();
					tmp.se(wd);
					data[j] = tmp;
					data[i].ae(wd);
				} else {
					ComplexField<N> wd = data[j];
					data[j] = data[i].s(wd);
					data[i].ae(wd);
				}
			}

			// a = 1 .. (dual-1) 
			for (int a = 1; a < dual; a++) {
				// trig recurrence for w-> exp(i theta) w 
				{
					if (MUTABLE) {
						N nn = n1;
						nn.se(s2);
						ComplexField<N> tmp = new ComplexField<N>(nn, s);
						w.me(tmp);
					} else {
						//w.me(c.fromPolar(n1, theta));
						w.me(new ComplexField<N>(n1.s(s2), s));
					}
				}
				for (int b = 0; b < n; b += 2 * dual) {
					int i = (b + a);
					int j = (b + a + (int) dual);

					ComplexField<N> z1 = data[j];

					if (MUTABLE) {
						ComplexField<N> wd = w.copy();
						wd.me(z1);

						ComplexField<N> tt = data[i].copy();
						tt.se(wd);
						data[j] = tt;
						data[i].ae(wd);
					} else {
						ComplexField<N> wd = w.m(z1);
						data[j] = data[i].s(wd);
						data[i].ae(wd);
					}
				}
			}*/


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

	protected void bitreverse(ComplexField<N>[] data) {
		/* This is the Goldrader bit-reversal algorithm */
		int n = data.length;
		int nm1 = n - 1;
		int i = 0;
		int j = 0;
		for (; i < nm1; i++) {

			// int ii = 2*i;
			// int ii = i << 1;

			// int jj = 2*j;
			// int jj = j << 1;

			// int k = n / 2 ;
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
	}
}
