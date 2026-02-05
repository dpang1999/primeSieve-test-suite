package generic;
public class GenSOR
{
	public static final double num_flops(int M, int N, int num_iterations)
	{
		double Md = (double) M;
		double Nd = (double) N;
		double num_iterD = (double) num_iterations;

		return (Md-1)*(Nd-1)*num_iterD*6.0;
	}

	public static final <R extends IField<R>> void execute(R omega, R G[][], int num_iterations)
	{
		int M = G.length;
		int N = G[0].length;

		R omega_over_four = omega.d(omega.coerce(4));
		R one_minus_omega = omega.coerce(1.0).s(omega);

		// update interior points
		//
		int Mm1 = M-1;
		int Nm1 = N-1; 
		for (int p=0; p<num_iterations; p++)
		{
			for (int i=1; i<Mm1; i++)
			{
				R[] Gi = G[i];
				R[] Gim1 = G[i-1];
				R[] Gip1 = G[i+1];
				for (int j=1; j<Nm1; j++)
					Gi[j] = omega_over_four.m(Gim1[j].a(Gip1[j]).a(Gi[j-1]).a(Gi[j+1])).a(one_minus_omega.m(Gi[j]));
			}
		}
	}

	// Main
	public static void main(String[] args)
	{
		int mode = 0;
		if (mode == 0) {
			int M = 10;
			int N = 10;
			int num_iterations = 100;

			DoubleField omega = new DoubleField(1.5);
			DoubleField G[][] = new DoubleField[M][N];
			for (int i=0; i<M; i++)
				for (int j=0; j<N; j++)
				{
					if (i == 0) {
						G[i][j] = new DoubleField(1.0); // Top edge is hot
					}
					else {
						G[i][j] = new DoubleField(0.0);
					}
				}
			System.out.println("Initial grid:");
			printMatrix(G);
					

			//System.out.println("Num flops: " + num_flops(M, N, num_iterations));
			execute(omega, G, num_iterations);

			System.out.println("Resulting grid:");
			printMatrix(G);
		}
		else {
			int M = 10;
			int N = 10;
			int num_iterations = 100;
			int prime = 40961;
			IntModP.setModulus(prime);
			IntModP omega = new IntModP(2);
			IntModP G[][] = new IntModP[M][N];
			for (int i=0; i<M; i++)
				for (int j=0; j<N; j++)
				{
					if (i == 0) {
						G[i][j] = new IntModP(1000); // Top edge is hot
					}
					else {
						G[i][j] = new IntModP(0);
					}
				}

			System.out.println("Initial grid:");
			printMatrix(G);
					

			//System.out.println("Num flops: " + num_flops(M, N, num_iterations));
			execute(omega, G, num_iterations);

			System.out.println("Resulting grid:");
			printMatrix(G);
		}
	}

	static <U> void printMatrix(U A[][]) {
		int M = A.length;
		int N = A[0].length;

		for (int i = 0; i < M; i++) {
			U Ai[] = A[i];
			for (int j = 0; j < N; j++)
				System.out.print(Ai[j].toString() + " ");
			System.out.println();
		}
		System.out.println();
	}

}
			
