package generic;
import java.util.Random;
public class GenMonteCarlo
{
	final static int SEED = 113;

	public static final double num_flops(int Num_samples)
	{
		// 3 flops in x^2+y^2 and 1 flop in random routine

		return ((double) Num_samples)* 4.0;

	}

	

	public static final <T extends IField<T>> T integrate(T t, int Num_samples)
	{

		Random R = new Random(SEED);

		int under_curve = 0;
		for (int count=0; count<Num_samples; count++)
		{
			T x= t.coerce(R.nextDouble());
			T y= t.coerce(R.nextDouble());

			if ( x.m(x).a(y.m(y)).coerce() <= 1.0)
				 under_curve ++;
			
		}

		return t.coerce(under_curve).d(t.coerce(Num_samples)).m(t.coerce(4.0));
	}

	public static void main(String[] args) {
		int mode = 0;
		if (mode == 0) {
			int Num_samples = 1000000;
			if (args.length > 0)
				Num_samples = Integer.parseInt(args[0]);
			DoubleField t = new DoubleField(0.0);
			DoubleField pi = integrate(t,Num_samples);
			System.out.println("Pi is approximately: " + pi);
			System.out.println("Num samples: " + Num_samples);
			System.out.println("Num flops: " + num_flops(Num_samples));
			System.out.println("RMS Error: " + Math.abs(Math.PI - pi.coerce()));
		}
		else {
			int Num_samples = 1000000;
			if (args.length > 0)
				Num_samples = Integer.parseInt(args[0]);
			int prime = 40961;
			IntModP.setModulus(prime);
			IntModP t = new IntModP(0);
			IntModP pi = integrate(t,Num_samples);
			System.out.println("Pi is approximately: " + pi);
			System.out.println("Num samples: " + Num_samples);
			System.out.println("Num flops: " + num_flops(Num_samples));
			System.out.println("RMS Error: " + Math.abs(Math.PI - pi.coerce()));
		}
	}

}
