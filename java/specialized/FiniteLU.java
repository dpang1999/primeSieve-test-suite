package specialized;
import helpers.LCG;
/**
    LU matrix factorization. (Based on TNT implementation.)
    Decomposes a matrix A  into a triangular lower triangular
    factor (L) and an upper triangular factor (U) such that
    A = L*U.  By convnetion, the main diagonal of L consists
    of 1's so that L and U can be stored compactly in
    a NxN matrix.


*/
public class FiniteLU 
{
    /**
        Returns a <em>copy</em> of the compact LU factorization.
        (useful mainly for debugging.)

        @return the compact LU factorization.  The U factor
        is stored in the upper triangular portion, and the L
        factor is stored in the lower triangular portion.
        The main diagonal of L consists (by convention) of
        ones, and is not explicitly stored.
    */
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

    protected static int[] new_copy(int x[])
    {
        int N = x.length;
        int T[] = new int[N];
        for (int i=0; i<N; i++)
            T[i] = x[i];
        return T;
    }


    protected static int[][] new_copy(int A[][])
    {
        int M = A.length;
        int N = A[0].length;

        int T[][] = new int[M][N];

        for (int i=0; i<M; i++)
        {
            int Ti[] = T[i];
            int Ai[] = A[i];
            for (int j=0; j<N; j++)
                Ti[j] = Ai[j];
        }

        return T;
    }




    protected static final void insert_copy(int B[][], int A[][])
    {
        int M = A.length;
        int N = A[0].length;

		int remainder = N & 3;		 // N mod 4;

        for (int i=0; i<M; i++)
        {
            int Bi[] = B[i];
            int Ai[] = A[i];
			for (int j=0; j<remainder; j++)
                Bi[j] = Ai[j];
            for (int j=remainder; j<N; j+=4)
			{
				Bi[j] = Ai[j];
				Bi[j+1] = Ai[j+1];
				Bi[j+2] = Ai[j+2];
				Bi[j+3] = Ai[j+3];
			}
		}
        
    }
    public int[][] getLU()
    {
        return new_copy(LU_);
    }

    /**
        Returns a <em>copy</em> of the pivot vector.

        @return the pivot vector used in obtaining the
        LU factorzation.  Subsequent solutions must
        permute the right-hand side by this vector.

    */
    public int[] getPivot()
    {
        return new_copy(pivot_);
    }
    
    /**
        Initalize LU factorization from matrix.

        @param A (in) the matrix to associate with this
                factorization.
    */
    public FiniteLU( int A[][], int modulus )
    {
        int M = A.length;
        int N = A[0].length;

        //if ( LU_ == null || LU_.length != M || LU_[0].length != N)
            LU_ = new int[M][N];

        insert_copy(LU_, A);

        //if (pivot_.length != M)
            pivot_ = new int[M];

        factor(LU_, pivot_, modulus);
    }

    /**
        Solve a linear system, with pre-computed factorization.

        @param b (in) the right-hand side.
        @return solution vector.
    */
    public int[] solve(int b[], int modulus)
    {
        int x[] = new_copy(b);

        solve(LU_, pivot_, x, modulus);
        return x;
    }
    

/**
    LU factorization (in place).

    @param A (in/out) On input, the matrix to be factored.
        On output, the compact LU factorization.

    @param pivit (out) The pivot vector records the
        reordering of the rows of A during factorization.
        
    @return 0, if OK, nozero value, othewise.
*/
public static int factor(int A[][],  int pivot[], int modulus)
{
    int N = A.length;
    int M = A[0].length;

    int minMN = Math.min(M,N);

    for (int j=0; j<minMN; j++)
    {
        // find pivot in column j and  test for singularity.

        int jp=j;
        
        double t = Math.abs(A[j][j]);
        for (int i=j+1; i<M; i++)
        {
            double ab = Math.abs(A[i][j]);
            if ( ab > t)
            {
                jp = i;
                t = ab;
            }
        }
        
        pivot[j] = jp;

        // jp now has the index of maximum element 
        // of column j, below the diagonal

        if ( A[jp][j] == 0 ) {
            System.out.println("Matrix is singular");                   
            return 1;       // factorization failed because of zero pivot
        }

        if (jp != j)
        {
            // swap rows j and jp
            int tA[] = A[j];
            A[j] = A[jp];
            A[jp] = tA;
        }

        if (j<M-1)                // compute elements j+1:M of jth column
        {
            // note A(j,j), was A(jp,p) previously which was
            // guarranteed not to be zero (Label #1)
            //
            int recp = modInverse(A[j][j], modulus);

            for (int k=j+1; k<M; k++)
                A[k][j] = A[k][j] * recp % modulus;
        }


        if (j < minMN-1)
        {
            // rank-1 update to trailing submatrix:   E = E - x*y;
            //
            // E is the region A(j+1:M, j+1:N)
            // x is the column vector A(j+1:M,j)
            // y is row vector A(j,j+1:N)


            for (int ii=j+1; ii<M; ii++)
            {
                int Aii[] = A[ii];
                int Aj[] = A[j];
                int AiiJ = Aii[j];
                for (int jj=j+1; jj<N; jj++)
                  Aii[jj] = (Aii[jj] - (AiiJ * Aj[jj])) % modulus;

            }
        }
    }

    return 0;
}


    /**
        Solve a linear system, using a prefactored matrix
            in LU form.


        @param LU (in) the factored matrix in LU form. 
        @param pivot (in) the pivot vector which lists
            the reordering used during the factorization
            stage.
        @param b    (in/out) On input, the right-hand side.
                    On output, the solution vector.
    */
    public static void solve(int LU[][], int pvt[], int b[], int modulus)
    {
        int M = LU.length;
        int N = LU[0].length;
        int ii=0;

        for (int i=0; i<M; i++)
        {
            int ip = pvt[i];
            int sum = b[ip];

            b[ip] = b[i];
            if (ii==0)
                for (int j=ii; j<i; j++)
                    sum = (sum - (LU[i][j] * b[j]) % modulus + modulus) % modulus;
            else 
                if (sum == 0)
                    ii = i;
            b[i] = sum;
        }

        for (int i=N-1; i>=0; i--)
        {
            int sum = b[i];
            for (int j=i+1; j<N; j++)
                sum = (sum - (LU[i][j] * b[j]) % modulus + modulus) % modulus;
            b[i] = (sum * modInverse(LU[i][i], modulus)) % modulus;
        }
    }               


    private int LU_[][];
    private int pivot_[];

    public static void main(String[] args) throws Exception {
        int N = 4;
        if (args.length > 0)
            N = Integer.parseInt(args[0]);
        int A[][] = new int[N][N];
        int b[] = new int[N];
        int modulus = (int)(Math.pow(2, 19) - 1);
        LCG rand = new LCG(12345, 1345, 16645, 1013904);
        for (int i=0; i<N; i++)
        {
            int row_sum = 0;
            for (int j=0; j<N; j++) {
                int val = rand.nextInt() % modulus;
                row_sum += val;
                A[i][j] = val;
            }
            A[i][i] = (row_sum + rand.nextInt() + 1) % modulus; // Ensure diagonal dominance
            b[i] = rand.nextInt() % modulus;
        }

        System.out.println("Java specialized finite field LU");
        System.out.println("Matrix size: " + N);
        //printMatrix(A);
        for (int i = 0; i < 10; i++) {
            int Acopy[][] = new_copy(A);
            int bcopy[] = new_copy(b);
            int pivot[] = new int[N];

            factor(Acopy, pivot, modulus);
            
            //only needed for debugging
            /* System.out.println("b: ");
            printVector(b); */
            solve(Acopy, pivot, bcopy, modulus); 
            /* System.out.println("Solution: ");
            printVector(b); */
            System.out.println("Iteration " + i + " completed");
        }
        

		System.exit(0);
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
    public static void printVector(double B[])
    {
        int N = B.length;

        for (int i=0; i<N; i++)
            System.out.print(B[i] + " ");
        System.out.println();
    }
}


