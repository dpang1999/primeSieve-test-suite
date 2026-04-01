package generic;

import java.lang.reflect.*;

import helpers.LCG;


/**
 * LU matrix factorization. (Based on TNT implementation.) Decomposes a matrix A
 * into a triangular lower triangular factor (L) and an upper triangular factor
 * (U) such that A = L*U. By convnetion, the main diagonal of L consists of 1's
 * so that L and U can be stored compactly in a NxN matrix.
 * 
 * 
 */
public class GenLU<R extends IField<R> & IOrdered<R> & IMath<R> & ICopiable<R>> {
	/**
	 * Returns a <em>copy</em> of the compact LU factorization. (useful mainly
	 * for debugging.)
	 * 
	 * @return the compact LU factorization. The U factor is stored in the upper
	 *         triangular portion, and the L factor is stored in the lower
	 *         triangular portion. The main diagonal of L consists (by
	 *         convention) of ones, and is not explicitly stored.
	 */


	@SuppressWarnings("unchecked")
	protected static <U extends IField<U>> U[] new_copy(U x[]) {
		int N = x.length;
		U T[] = (U[]) Array.newInstance(x[0].getClass(), N);
		for (int i = 0; i < N; i++)
			T[i] = x[i];
		return T;
	}

	@SuppressWarnings("unchecked")
	protected static <U extends IField<U>> U[][] new_copy(U A[][]) {
		int M = A.length;
		int N = A[0].length;
		int dims[] = new int[2];
		dims[0] = M;
		dims[1] = N;

		U T[][] = (U[][]) Array.newInstance(A[0][0].getClass(), dims);// double[M][N]

		for (int i = 0; i < M; i++) {
			U Ti[] = T[i];
			U Ai[] = A[i];
			for (int j = 0; j < N; j++)
				Ti[j] = Ai[j];
		}

		return T;
	}

	public static int[] new_copy(int x[]) {
		int N = x.length;
		int T[] = new int[N];
		for (int i = 0; i < N; i++)
			T[i] = x[i];
		return T;
	}

	protected static final <U extends IField<U>> void insert_copy(U B[][],
			U A[][]) {
		int M = A.length;
		int N = A[0].length;

		int remainder = N & 3; // N mod 4;

		for (int i = 0; i < M; i++) {
			U Bi[] = B[i];
			U Ai[] = A[i];
			for (int j = 0; j < remainder; j++)
				Bi[j] = Ai[j];
			for (int j = remainder; j < N; j += 4) {
				Bi[j] = Ai[j];
				Bi[j + 1] = Ai[j + 1];
				Bi[j + 2] = Ai[j + 2];
				Bi[j + 3] = Ai[j + 3];
			}
		}
	}

	public R[][] getLU() {
		return new_copy(LU_);
	}

	/**
	 * Returns a <em>copy</em> of the pivot vector.
	 * 
	 * @return the pivot vector used in obtaining the LU factorzation.
	 *         Subsequent solutions must permute the right-hand side by this
	 *         vector.
	 * 
	 */
	public int[] getPivot() {
		return new_copy(pivot_);
	}

	/**
	 * Initalize LU factorization from matrix.
	 * 
	 * @param A
	 *            (in) the matrix to associate with this factorization.
	 */
	@SuppressWarnings("unchecked")
	public GenLU(R A[][]) {
		int M = A.length;
		int N = A[0].length;
		int dims[] = new int[2];
		dims[0] = M;
		dims[1] = N;

		// if ( LU_ == null || LU_.length != M || LU_[0].length != N)
		LU_ = (R[][]) Array.newInstance(A[0][0].getClass(), dims);// double[M][N]

		insert_copy(LU_, A);

		// if (pivot_.length != M)
		pivot_ = new int[M];

		factor(LU_, pivot_);
	}

	/**
	 * Solve a linear system, with pre-computed factorization.
	 * 
	 * @param b
	 *            (in) the right-hand side.
	 * @return solution vector.
	 */
	public R[] solve(R b[]) {
		R x[] = new_copy(b);

		solve(LU_, pivot_, x);
		return x;
	}

	/**
	 * LU factorization (in place).
	 * 
	 * @param A
	 *            (in/out) On input, the matrix to be factored. On output, the
	 *            compact LU factorization.
	 * 
	 * @param pivit
	 *            (out) The pivot vector records the reordering of the rows of A
	 *            during factorization.
	 * 
	 * @return 0, if OK, nozero value, othewise.
	 */
	public static <U extends IField<U> & IOrdered<U> & IMath<U> & ICopiable<U>> int factor(U A[][], int pivot[]) {
		int N = A.length;
		int M = A[0].length;

		int minMN = Math.min(M, N);

		for (int j = 0; j < minMN; j++) {
			// find pivot in column j and test for singularity.

			int jp = j;

			U t = A[j][j].copy();
			t = t.abs();

			for (int i = j + 1; i < M; i++) {
				U ab = A[i][j].copy();
				ab = ab.abs();
				if (ab.coerce() > t.coerce()) {
					jp = i;
					t = ab;
				}
			}

			pivot[j] = jp;

			// jp now has the index of maximum element
			// of column j, below the diagonal

			if (A[jp][j].coerce() == 0) 
			{
				System.out.println("Matrix is singular");
				return 1; // factorization failed because of zero pivot
			}

			if (jp != j) {
				// swap rows j and jp
				U tA[] = A[j];
				A[j] = A[jp];
				A[jp] = tA;
			}

			if (j < M - 1) { // compute elements j+1:M of jth column
				// note A(j,j), was A(jp,p) previously which was
				// guarranteed not to be zero (Label #1)
				//
				U recp = A[j][j].one().d(A[j][j]); 

				for (int k = j + 1; k < M; k++)
					A[k][j] = A[k][j].m(recp);
			}

			if (j < minMN - 1) {
				// rank-1 update to trailing submatrix: E = E - x*y;
				//
				// E is the region A(j+1:M, j+1:N)
				// x is the column vector A(j+1:M,j)
				// y is row vector A(j,j+1:N)

				for (int ii = j + 1; ii < M; ii++) {
					U Aii[] = A[ii];
					U Aj[] = A[j];
					U AiiJ = Aii[j];
					for (int jj = j + 1; jj < N; jj++){
						U temp = AiiJ.m(Aj[jj]);
						Aii[jj] = Aii[jj].s(temp);
					}
				}
			}
		}

		return 0;
	}

	/**
	 * Solve a linear system, using a prefactored matrix in LU form.
	 * 
	 * 
	 * @param LU
	 *            (in) the factored matrix in LU form.
	 * @param pivot
	 *            (in) the pivot vector which lists the reordering used during
	 *            the factorization stage.
	 * @param b
	 *            (in/out) On input, the right-hand side. On output, the
	 *            solution vector.
	 */
	public static <U extends IField<U>> void solve(U LU[][], int pvt[], U b[]) {
		int M = LU.length;
		int N = LU[0].length;
		int ii = 0;

		for (int i = 0; i < M; i++) {
			int ip = pvt[i];
			U sum = b[ip];

			b[ip] = b[i];
			if (ii == 0)
				for (int j = ii; j < i; j++){
					U temp = LU[i][j].m(b[j]);
					sum = sum.s(temp);
				}
			else if (sum.coerce() == 0.0)
				ii = i;
			b[i] = sum;
		}

		for (int i = N - 1; i >= 0; i--) {
			U sum = b[i];
			for (int j = i + 1; j < N; j++) {
				U temp = LU[i][j].m(b[j]);
				sum = sum.s(temp);
			}
				
			b[i] = sum.d(LU[i][i]);
		}
	}

	private R LU_[][];

	private int pivot_[];

	public static void main(String[] args) throws Exception {
		int N = 4;
		if (args.length > 0)
			N = Integer.parseInt(args[0]);
		int mode = 0;
		if (args.length > 1)
			mode = Integer.parseInt(args[1]);
		LCG rand = new LCG(12345, 1345, 16645, 1013904);
		if (mode == 0) {
			System.out.println("Java Generic DoubleField LU");
			System.out.println("Matrix size: " + N);
			DoubleField A[][] = new DoubleField[N][N];
			DoubleField b[] = new DoubleField[N];
			
			

			for (int i = 0; i < N; i++) {
				double row_sum = 0;
				for (int j = 0; j < N; j++)
				{
					double val = rand.nextDouble() * 1000;
					row_sum += val;
					A[i][j] = new DoubleField(val);
				}
				A[i][i] = new DoubleField(row_sum + rand.nextDouble()*1000 + 1.0); // Ensure diagonal dominance
				b[i] = new DoubleField(rand.nextDouble()*1000);
			}

			for (int i = 0; i < 10; i++) {
				//printMatrix(A);
				DoubleField ACopy[][] = new_copy(A);
				int pivot[] = new int[N];
				
				factor(ACopy, pivot);
				
				//System.out.println("b: ");
				//printVector(b);
				DoubleField BCopy[] = new_copy(b);
				solve(ACopy, pivot, BCopy);
				//System.out.println("Solution: ");
				//printVector(b);

				//DoubleField product[] = multiplyMatrices(A, BCopy);
				//printVector(product);

				//System.out.println("RMS Difference: " + RMSDiff(b, product));
				System.out.println("Iteration " + i + " completed");
				}
		}
		else {
			System.out.println("Java Generic IntModP LU");
			System.out.println("Matrix size: " + N);
			IntModP A[][] = new IntModP[N][N];
			IntModP b[] = new IntModP[N];
			int pivot[] = new int[N];
			long p = (long)(Math.pow(2, 13) - 1);
			IntModP.setModulus(p);

			for (int i = 0; i < N; i++) {
				int row_sum = 0;
				for (int j = 0; j < N; j++) {
					int value = rand.nextInt();
					A[i][j] = new IntModP(value);
					row_sum += value;
				}
				A[i][i]= new IntModP(row_sum + rand.nextInt() + 1); // Ensure diagonal dominance					
				b[i] = new IntModP(rand.nextInt());
			}

			//printMatrix(A);
				for (int i = 0; i < 10; i++) {
				IntModP ACopy[][] = new_copy(A);

				factor(ACopy, pivot);
				//System.out.println("b: ");
				//printVector(b);
				IntModP BCopy[] = new_copy(b);
				solve(ACopy, pivot, BCopy);
				/* System.out.println("Solution: ");
				printVector(b);

				IntModP product[] = multiplyMatrices(ACopy, b);
				printVector(product);

				System.out.println("RMS Difference: " + RMSDiff(BCopy, product)); */
				System.out.println("Iteration " + i + " completed");
			}
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

	static <U> void printVector(U B[]) {
		int M = B.length;

		for (int i = 0; i < M; i++)
			System.out.print(B[i].toString() + " ");
		System.out.println();
		System.out.println();
	}

	static <U extends IField<U>> U[] multiplyMatrices(U A[][], U B[]) {
		int M = A.length;
		int N = A[0].length;
		int P = B.length;

		if (N != P) {
			throw new IllegalArgumentException("Incompatible matrix dimensions");
		}

		// Multiply an NxN matrix by an Nx1 vector
		U C[] = (U[]) Array.newInstance(A[0][0].getClass(), M);
		for (int i = 0; i < M; i++) {
			C[i] = B[0].zero();
			for (int j = 0; j < N; j++) {
				C[i] = C[i].a(A[i][j].m(B[j]));
			}
		}
		return C;
	}

	static <U extends IField<U>> double RMSDiff(U A[], U B[]) {
		if (A.length != B.length) {
			throw new IllegalArgumentException("Incompatible vector dimensions");
		}

		double sum = 0;
		for (int i = 0; i < A.length; i++) {
			double diff = A[i].coerce() - B[i].coerce();
			sum += diff * diff;
		}
		return Math.sqrt(sum / A.length);
	}
}

