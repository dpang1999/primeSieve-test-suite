package generic;
import java.util.*;
import java.lang.reflect.*;
import helpers.LCG;
import helpers.FindPrime;

/**
 * A generic FFT computation for complex fields or finite fields 
 * where the basefield of the complex fields can be any generic field
 * 
 * @author Daniel Pang daniel.pang@uwaterloo.ca,
 * @author Derived from Laurentiu Dragan (ldragan@scl.csd.uwo.ca) 2005 Thesis
 */



public class GenFFT<N extends IField<N> & IOrdered<N> & ICopiable<N> & IPrimitiveRoots<N> & IMath<N>> {

	private N c;

	final static boolean DEBUG = false;
	final static boolean MUTABLE = false;

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
		N norm = c.coerce(nd);
		for (int i = 0; i < nd; i++) {
			//System.out.print(data[i]);
			data[i].de(norm);
			//System.out.println(" => " + data[i]);
		}
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

		// Transform & invert
		transform(data);

		System.out.println("After transform:" + Arrays.toString(data));

		inverse(data);

		System.out.println("After inverse:" + Arrays.toString(data));
		// Compute RMS difference. 
		double diff = 0.0;
		for (int i = 0; i < nd; i++) {
			N d = data[i].copy();
			N difference = d.s(copy.get(i));
			difference.abs();
			diff += difference.coerce()*difference.coerce();
		}

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
		int type = 2;
		if (args.length == 0 && type == 0) {
			int n = 4;
			ComplexField<DoubleField>[] data = fft.makeRandom(n);

			System.out.println(Arrays.toString(data));
			System.out.println("n=" + n + " => RMS Error="+ fft.test(data));


			//double[] in1 = {38.0,  0.0, 44.0, 87.0,  6.0, 45.0, 22.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			//double[] in2 = {80.0, 18.0, 62.0, 90.0, 17.0, 96.0, 27.0, 97.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			//double[] in1 = {11400, 28374, 23152, 9576, 29511, 20787, 13067, 14015, 0, 0, 0, 0, 0, 0, 0, 0};
			//double[] in2 = {30268, 20788, 8033, 15446, 26275, 11619, 2494, 7016, 0, 0, 0, 0, 0, 0, 0, 0};
			double[] in1 = {33243586, 638827078, 767661659, 778933286, 790244973, 910208076, 425757125,
				478004096, 153380495, 205851834, 668901196, 15731080, 899763115, 551605421,
				181279081, 600279047, 711828654, 483031418, 737709105, 20544909, 609397212,
				201989947, 215952988, 206613081, 471852626, 889775274, 992608567, 947438771,
				969970961, 676943009, 934992634, 922939225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
			double[] in2 = {194132110, 219972873, 66644114, 902841100, 565039275, 540721923, 810650854,
				702680360, 147944788, 859947137, 59055854, 288190067, 537655879, 836782561,
				308822170, 315498953, 417177801, 640439652, 198304612, 525827778, 115633328,
				285831984, 136721026, 203065689, 884961191, 222965182, 735241234, 746745227,
				667772468, 739110962, 610860398, 965331182, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

			
			ComplexField<DoubleField>[] data1 = (ComplexField<DoubleField>[]) Array.newInstance(c.getClass(), in1.length / 2);
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
			IntModP basicNum = new IntModP(1);
			GenFFT basicFft = new GenFFT(basicNum);
			IntModP[] basicData = new IntModP[basic.length];
			for(int i = 0; i < basic.length; i ++) {
				basicData[i] = basicNum.coerce(basic[i]);
			}
			System.out.println("n=" + basic.length + " => RMS Error="+ basicFft.test(basicData));


			System.out.println("Basic Data: " + Arrays.toString(basicData));
			basicFft.transform(basicData);
			System.out.println("Transformed Basic Data: " + Arrays.toString(basicData));
			basicFft.inverse(basicData);
			System.out.println("Inverse Basic Data: " + Arrays.toString(basicData));

			
			int[] in1 = {
				38,  0, 44, 87,  6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0
			};
			int[] in2 = {
				80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0
			};
			long[] out = {3040, 684, 5876, 11172,  5420, 16710, 12546, 20555,16730, 15704, 21665,  5490, 13887,  4645,  9021,0};
			long prime = 40961;
			IntModP.setModulus(prime);
			
			
			/*
			int[] in1 = {11400, 28374, 23152, 9576, 29511, 20787, 13067, 14015, 0, 0, 0, 0, 0, 0, 0, 0
			};
			int[] in2 = {30268, 20788, 8033, 15446, 26275, 11619, 2494, 7016, 0, 0, 0, 0, 0, 0, 0, 0
			};
			long[] out = { 345055200L, 1095807432L, 1382179648L, 1175142886L, 2016084656L, 2555168834L, 2179032777L, 1990011337L, 1860865174L, 1389799087L, 942120918L, 778961552L, 341270975L, 126631482L, 98329240L, 0L };
			long prime = 3221225473L;
			*/
			
			
			/*
			int[] in1 = {33243586, 638827078, 767661659, 778933286, 790244973, 910208076, 425757125,
				478004096, 153380495, 205851834, 668901196, 15731080, 899763115, 551605421,
				181279081, 600279047, 711828654, 483031418, 737709105, 20544909, 609397212,
				201989947, 215952988, 206613081, 471852626, 889775274, 992608567, 947438771,
				969970961, 676943009, 934992634, 922939225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
			int[] in2 = {194132110, 219972873, 66644114, 902841100, 565039275, 540721923, 810650854,
				702680360, 147944788, 859947137, 59055854, 288190067, 537655879, 836782561,
				308822170, 315498953, 417177801, 640439652, 198304612, 525827778, 115633328,
				285831984, 136721026, 203065689, 884961191, 222965182, 735241234, 746745227,
				667772468, 739110962, 610860398, 965331182, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
			long prime = 4179340454199820289L;
			*/
					
				
			
			
			IntModP finiteNum = new IntModP(1);
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
		else {
			LCG rand = new LCG(12345, 1345, 65, 17);
			/*int[] randomNumbers = new int[10];
			double[] randomDoubles = new double[10];
			for (int i = 0; i < randomNumbers.length; ic ++) {
				randomNumbers[i] = rand.nextInt();
				randomDoubles[i] = rand.nextDouble();
			}
			System.out.println("Random Numbers: " + Arrays.toString(randomNumbers));
			System.out.println("Random Doubles: " + Arrays.toString(randomDoubles));*/
		
			int n = Integer.parseInt(args[0]);
			int fieldType = Integer.parseInt(args[1]); // 0 for finite field, 1 for complex field
			if (fieldType == 0) {
				IntModP[] data = new IntModP[n];
				int prime = FindPrime.findPrimeCongruentOneModN(n);
				IntModP.setModulus(prime);
				for (int i = 0; i < n; i++) {
					data[i] = new IntModP(rand.nextInt() % prime);
				}
				GenFFT<IntModP> finiteFft = new GenFFT<IntModP>(new IntModP(0));
				System.out.println("Generic Java FFT Tests");
				System.out.println("Java Generics, Finite Field, n=" + n);
				for (int i = 0; i < 10; i++) {
					finiteFft.transform_internal(data, -1);
					finiteFft.transform_internal(data, 1);
					System.out.println("Loop " + i + " done.");	
				}

				
			}
			else {
				ComplexField<DoubleField>[] data = (ComplexField<DoubleField>[]) Array.newInstance(c.getClass(), n);
				for (int i = 0; i < n; i++) {
					data[i] = new ComplexField<DoubleField>(
						new DoubleField(rand.nextDouble()),
						new DoubleField(rand.nextDouble()));
				}
				ComplexField<DoubleField> c1 = new ComplexField<DoubleField>(new DoubleField(0),
					new DoubleField(0));
				GenFFT<ComplexField<DoubleField>> fft1 = new GenFFT<ComplexField<DoubleField>>(c1);
				System.out.println("Generic Java FFT Tests");
				System.out.println("Java Generics, Complex Field, n=" + n);				
				for (int i = 0; i<= 10; i++) {
					fft1.transform_internal(data, -1);
					fft1.transform_internal(data, 1);
					System.out.println("Loop " + i + " done.");
				}
			}
		
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


	protected void transform_internal(N[] data, int direction) {
		if (data.length == 0)
			return;
		int n = data.length;
		if (n == 1)
			return; // Identity operation
		int logn = log2(n);

		// bit reverse
		bitreverse(data);
		//System.out.println("After Bit-Reversal (Java): " + Arrays.toString(data));

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
			}
			//System.out.println("After Stage " + bit + " (Java): " + Arrays.toString(data));
		}


	}

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
