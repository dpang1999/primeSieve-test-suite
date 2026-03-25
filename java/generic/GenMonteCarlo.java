package generic;
import helpers.LCG;
public class GenMonteCarlo
{

	public static final <T extends IField<T> & IOrdered<T>> double integrate(T t, int Num_samples)
	{
		LCG rand = new LCG(12345, 1345, 16645, 1013904);
		int under_curve = 0;
		for (int count=0; count<Num_samples; count++)
		{
			T x= t.coerce(rand.nextDouble());
			T y= t.coerce(rand.nextDouble());

			if ( x.m(x).a(y.m(y)).le(t.one()))
				 under_curve ++;
			
		}

		return (double) under_curve / (double) Num_samples;
	}

	public static void main(String[] args) {
		int mode = 0;
		if (args.length > 0)
			mode = Integer.parseInt(args[0]);
		if (mode == 0) {
			int Num_samples = 1000000;
			if (args.length > 1)
				Num_samples = Integer.parseInt(args[1]);
			DoubleField t = new DoubleField(0.0);
			double pi = integrate(t,Num_samples);
			System.out.println("Java generic doublefield montecarlo");
			System.out.println("Pi is approximately: " + pi);
			System.out.println("Num samples: " + Num_samples);
			System.out.println("RMS Error: " + Math.abs(Math.PI - pi));
		}
		else if (mode == 1) {
			int Num_samples = 1000000;
			if (args.length > 1)
				Num_samples = Integer.parseInt(args[1]);
			SingleField t = new SingleField(0.0f);
			double pi = integrate(t,Num_samples);
			System.out.println("Java generic singlefield montecarlo");
			System.out.println("Pi is approximately: " + pi);
			System.out.println("Num samples: " + Num_samples);
			System.out.println("RMS Error: " + Math.abs(Math.PI - pi));
		}
		else {
			int Num_samples = 1000000;
			if (args.length > 1)
				Num_samples = Integer.parseInt(args[1]);
			int prime = 40961;
			IntModP.setModulus(prime);
			IntModP t = new IntModP(0);
			double pi = integrate(t,Num_samples);
			System.out.println("Java generic intmodp montecarlo");
			System.out.println("Pi is approximately: " + pi);
			System.out.println("Num samples: " + Num_samples);
			System.out.println("RMS Error: " + Math.abs(Math.PI - pi));
		}
	}

}
