package specialized;

public class SOR
{
	public static final double num_flops(int M, int N, int num_iterations)
	{
		double Md = (double) M;
		double Nd = (double) N;
		double num_iterD = (double) num_iterations;

		return (Md-1)*(Nd-1)*num_iterD*6.0;
	}

	public static final void execute(double omega, double G[][], int 
			num_iterations)
	{
		int M = G.length;
		int N = G[0].length;

		double omega_over_four = omega * 0.25;
		double one_minus_omega = 1.0 - omega;

		// update interior points
		//
		int Mm1 = M-1;
		int Nm1 = N-1; 
		for (int p=0; p<num_iterations; p++)
		{
			for (int i=1; i<Mm1; i++)
			{
				double[] Gi = G[i];
				double[] Gim1 = G[i-1];
				double[] Gip1 = G[i+1];
				for (int j=1; j<Nm1; j++)
					Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j-1] 
								+ Gi[j+1]) + one_minus_omega * Gi[j];
			}
		}
	}

			
    public static void main(String[] args)
    {
        int M = 10;
        int N = 10;
        int num_iterations = 1;

        double omega = 1.5;
        double[][] G = new double[M][N];

        // Set boundary conditions
        for (int i = 0; i < M; i++) {
            G[i][0] = 0;         // Left edge
            G[i][N-1] = 0;       // Right edge
        }
        for (int j = 0; j < N; j++) {
            G[0][j] = 100;       // Top edge (hot)
            G[M-1][j] = 0;       // Bottom edge (cold)
        }

        printMatrix(G);
        SOR.execute(omega, G, num_iterations);


        System.out.println("\nSteady-state temperature distribution:");
        printMatrix(G);
        
    }

    public static void printMatrix(double A[][])
    {
        int M = A.length;
        int N = A[0].length;

        for (int i=0; i<M; i++)
        {
            for (int j=0; j<N; j++)
                System.out.print(A[i][j] + " ");
            System.out.println();
        }
    }
}