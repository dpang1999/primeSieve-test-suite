package generic;

import java.lang.reflect.*;


/**
 * LU matrix factorization. (Based on TNT implementation.) Decomposes a matrix A
 * into a triangular lower triangular factor (L) and an upper triangular factor
 * (U) such that A = L*U. By convnetion, the main diagonal of L consists of 1's
 * so that L and U can be stored compactly in a NxN matrix.
 * 
 * 
 */
public class GenLU<R extends IRing<R> & IInvertible<R>> {
	/**
	 * Returns a <em>copy</em> of the compact LU factorization. (useful mainly
	 * for debugging.)
	 * 
	 * @return the compact LU factorization. The U factor is stored in the upper
	 *         triangular portion, and the L factor is stored in the lower
	 *         triangular portion. The main diagonal of L consists (by
	 *         convention) of ones, and is not explicitly stored.
	 */

	public static final double num_flops(int N) {
		// rougly 2/3*N^3

		double Nd = (double) N;
		return (2.0 * Nd * Nd * Nd / 3.0);
	}

	@SuppressWarnings("unchecked")
	protected static <U extends IRing<U>> U[] new_copy(U x[]) {
		int N = x.length;
		U T[] = (U[]) Array.newInstance(x[0].getClass(), N);
		for (int i = 0; i < N; i++)
			T[i] = x[i];
		return T;
	}

	@SuppressWarnings("unchecked")
	protected static <U extends IRing<U>> U[][] new_copy(U A[][]) {
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

	protected static final <U extends IRing<U>> void insert_copy(U B[][],
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
	public static <U extends IRing<U> & IInvertible<U>> int factor(U A[][], int pivot[]) {
		int N = A.length;
		int M = A[0].length;

		int minMN = Math.min(M, N);

		for (int j = 0; j < minMN; j++) {
			// find pivot in column j and test for singularity.

			int jp = j;

			U t = A[j][j].coerce(Math.abs(A[j][j].coerce()));
			for (int i = j + 1; i < M; i++) {
				U ab = A[i][j].coerce(Math.abs(A[i][j].coerce()));
				if (ab.coerce() > t.coerce()) {
					jp = i;
					t = ab;
				}
			}

			pivot[j] = jp;

			// jp now has the index of maximum element
			// of column j, below the diagonal

			if (A[jp][j].coerce() == 0)
				return 1; // factorization failed because of zero pivot

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
				U recp = A[j][j].coerce(1.0).m(A[j][j].invert());

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
					for (int jj = j + 1; jj < N; jj++)
						Aii[jj] = Aii[jj].s(AiiJ.m(Aj[jj]));
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
	public static <U extends IRing<U> & IInvertible<U>> void solve(U LU[][], int pvt[], U b[]) {
		int M = LU.length;
		int N = LU[0].length;
		int ii = 0;

		for (int i = 0; i < M; i++) {
			int ip = pvt[i];
			U sum = b[ip];

			b[ip] = b[i];
			if (ii == 0)
				for (int j = ii; j < i; j++)
					sum = sum.s(LU[i][j].m(b[j]));
			else if (sum.coerce() == 0.0)
				ii = i;
			b[i] = sum;
		}

		for (int i = N - 1; i >= 0; i--) {
			U sum = b[i];
			for (int j = i + 1; j < N; j++)
				sum = sum.s(LU[i][j].m(b[j]));
			b[i] = sum.m(LU[i][i].invert());
		}
	}

	private R LU_[][];

	private int pivot_[];

	public static void main(String[] args) throws Exception {
		int N = 4;
		if (args.length > 0)
			N = Integer.parseInt(args[0]);

		DoubleRing A[][] = new DoubleRing[N][N];
		DoubleRing b[] = new DoubleRing[N];
		int pivot[] = new int[N];
		

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				A[i][j] = new DoubleRing(Math.random()*1000);
			b[i] = new DoubleRing(Math.random()*1000);
		}
		printMatrix(A);

		factor(A, pivot);
		
		System.out.println("b: ");
		printVector(b);
		//solve(A, pivot, b); //only needed for debugging
		System.out.println("Solution: ");
		printVector(b);
		System.exit(0);
	}

	static void printMatrix(DoubleRing A[][]) {
		int M = A.length;
		int N = A[0].length;

		for (int i = 0; i < M; i++) {
			DoubleRing Ai[] = A[i];
			for (int j = 0; j < N; j++)
				System.out.print(Ai[j].toString() + " ");
			System.out.println();
		}
		System.out.println();
	}
	static void printVector(DoubleRing B[]) {
		int M = B.length;

		for (int i = 0; i < M; i++)
			System.out.print(B[i].toString() + " ");
		System.out.println();
		System.out.println();
	}
}

