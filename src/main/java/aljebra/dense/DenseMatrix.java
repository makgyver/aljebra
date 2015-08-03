/*******************************************************************************
 * Copyright (C) 2015 Mirko Polato
 *
 * This file is part of aljebra.
 *
 * aljebra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * aljebra is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with aljebra. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package aljebra.dense;

import java.io.Serializable;
import java.util.ArrayList;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import aljebra.sparse.SparseMatrix;
import aljebra.utils.Misc;
import aljebra.utils.RNG;

/**
 * <p>Dense implementation of a matrix.</p>
 * 
 * <strong>This implementation is a re-elaboration of the one
 * provided by Guibing Guo in his <a href="http://www.librec.net">librec</a>.</strong>
 * 
 * @author Mirko Polato
 *
 */
public class DenseMatrix implements Cloneable, Serializable {

	private static final long serialVersionUID = -989874464674003990L;
	
	//
	// STATIC SECTION
	//
	
	/**
	 * Constructs a new zero matrix with dimension 
	 * {@code rows} x {@code cols}.
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 * @return a new dense null matrix
	 */
	public static DenseMatrix zero(int rows, int cols) {
		return new DenseMatrix(rows, cols);
	}
	
	/**
	 * Constructs a new matrix with all entries set to one,
	 * with dimension {@code rows} x {@code cols}.
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 * @return a new one matrix
	 */
	public static DenseMatrix one(int rows, int cols) {
		DenseMatrix one = new DenseMatrix(rows, cols);
		one.setAll(1.0);
		return one;
	}
	
	/**
	 * Constructs a new matrix with random entries in the range
	 * of real [0,1], with dimension {@code rows} x {@code cols}.
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 * @return a new random matrix
	 */
	public static DenseMatrix rand(int rows, int cols) {
		DenseMatrix vec = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				vec.set(i, j, RNG.uniformDbl());
			}
		}
		return vec;
	}
	
	/**
	 * Constructs an identity matrix with dimension 
	 * {@code n} x {@code n}.
	 * 
	 * @param n	 number of rows and column
	 * @return a new identity matrix
	 */
	public static DenseMatrix identity(int n) {
		DenseMatrix vec = new DenseMatrix(n, n);
		for (int i = 0; i < n; ++i) {
			vec.set(i, i, 1.0);
		}
		return vec;
	}
	
	/**
	 * Performs the dot product between the a row of the {@code M} matrix
	 * and a column of the {@code N} matrix. The dimensions of the matrices
	 * should be <tt>n x m</tt> and <tt>m x p</tt>.
	 * 
	 * @param M		the first matrix
	 * @param mrow	the row of the first matrix
	 * @param N		the second matrix
	 * @param ncol	the column of the second matrix
	 * @return the M-row by N-column dot product 
	 */
	public static double dot(DenseMatrix M, int mrow, DenseMatrix N, int ncol) {
		assert M.cols == N.rows;

		double result = 0;
		for (int j = 0; j < M.cols; ++j) {
			result += M.get(mrow, j) * N.get(j, ncol);
		}
		return result;
	}
	
	/**
	 * Constructs a new diagonal matrix, where the values
	 * in the diagonal are taken from the given vector.
	 * 
	 * @param vec	the diagonal vector
	 * @return a new diagonal matrix
	 */
	public static DenseMatrix diagonal(DenseVector vec) {
		DenseMatrix diag = new DenseMatrix(vec.size(), vec.size());
		for (int i = 0; i < diag.rows; ++i) {
			diag.data[i][i] = vec.get(i);
		}
		return diag;
	}
	
	/**
	 * Fast implementation of the dot product <<tt>A.T, A</tt>>.
	 * 
	 * @param matrix	the matrix
	 * @return the dot product
	 */
	public static DenseMatrix AtA(DenseMatrix matrix) {
		DenseMatrix res = new DenseMatrix(matrix.rows, matrix.rows);

		for (int i = 0; i < matrix.rows; i++) {
			for (int k = 0; k < matrix.rows; k++) {
				double val = 0;
				for (int j = 0; j < matrix.cols; j++) {
					val += matrix.data[i][j] * matrix.data[k][j];
				}
				res.data[i][k] = val;
			}
		}
		return res;
	}
	
	//
	// END STATIC SECTION
	//
	
	// Number of rows
	protected final int rows;
	
	// Number of columns
	protected final int cols;
	
	// The actual data
	protected double[][] data;
	
	/**
	 * Constructs a new dense matrix with dimension 
	 * {@code rows} x {@code cols}.
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 */
	public DenseMatrix(int rows, int cols) {
		this.rows = rows;
		this.cols = cols;
		this.data = new double[rows][cols]; 
	}
	
	/**
	 * Constructs a new dense matrix using the given
	 * bidimensional {@code array} as data. The
	 * dimension of the matrix is inferred from the array.
	 * 
	 * @param array	 the data array 
	 */
	public DenseMatrix(double[][] array) {
		this.rows = array.length;
		this.cols = array[0].length;
		this.data = new double[rows][cols];
		
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				this.data[i][j] = array[i][j];
			}
		}
	}
	
	/**
	 * Constructs a new dense matrix using the given
	 * bidimensional {@code array} as data. The
	 * dimension of the matrix will be {@code rows} x {@code cols}.
	 * 
	 * @param array	 the data array
	 * @param rows	 number of rows
	 * @param cols	 number of columns
	 */
	public DenseMatrix(double[][] array, int rows, int cols) {
		this.rows = rows;
		this.cols = cols;
		this.data = array;
	}
	
	/**
	 * Copy constructor.
	 * 
	 * @param mat	the dense matrix to copy
	 */
	public DenseMatrix(DenseMatrix mat) {
		this(mat.data);
	}
	
	@Override
	public DenseMatrix clone() {
		return new DenseMatrix(this);
	}

	/**
	 * Returns the matrix entry in position ({@code row}, {@code col}).
	 * 
	 * @param row	the row of the entry
	 * @param col	the column of the entry
	 * @return	the entry in position {@code row}, {@code col}
	 */
	public double get(int row, int col) {
		assert(row >= 0 && row < rows && col >= 0 && col < cols);
		return data[row][col];
	}
	
	/**
	 * Sets the entry in position ({@code row}, {@code col}) to the given {@code value}.
	 * 
	 * @param row	 the row of the entry
	 * @param col	 the column of the entry
	 * @param value	 the new value for the entry
	 */
	public void set(int row, int col, double value) {
		assert(row >= 0 && row < rows && col >= 0 && col < cols);
		data[row][col] = value;
	}
	
	/**
	 * Sets all the entries to the given {@code value}.
	 * 
	 * @param value	 the new value for the entries
	 */
	public void setAll(double value) {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				data[i][j] = value;
			}
		}
	}
	
	/**
	 * Returns the number of rows of the matrix.
	 * 
	 * @return the number of rows of the matrix
	 */
	public int rows() {
		return rows;
	}
	
	/**
	 * Returns the number of columns of the matrix.
	 * 
	 * @return the number of columns of the matrix
	 */
	public int cols() {
		return cols;
	}
	
	/**
	 * Returns true if the matrix is square.
	 * 
	 * @return whether the matrix is square or not
	 */
	public boolean isSquare() {
		return rows == cols;
	}
	
	/**
	 * Returns the {@code i}-th row as a dense vector.
	 * 
	 * @param i	 the row's index
	 * @return the {@code i}-th row as a dense vector
	 */
	public DenseVector row(int i) {
		assert(i >= 0 && i < rows);
		return new DenseVector(data[i]);
	}
	
	/**
	 * Returns the {@code j}-th column as a dense vector.
	 * 
	 * @param j  the column's index
	 * @return the {@code j}-th column as a dense vector
	 */
	public DenseVector col(int j) {
		assert(j >= 0 && j < cols);
		
		DenseVector result = new DenseVector(rows);
		for (int i = 0; i < rows; ++i) {
			result.set(i, data[i][j]);
		}
		return result;
	}
	
	/**
	 * Returns the main diagonal of the matrix as a 
	 * dense vector.
	 * 
	 * @return the main diagonal of the matrix as a dense vector
	 */
	public DenseVector diag() {
		assert(isSquare());
		
		DenseVector result = new DenseVector(rows);
		for (int i = 0; i < rows; ++i) {
			result.set(i, data[i][i]);
		}
		return result;
	}
	
	/**
	 * Changes a row with the given vector.
	 * 
	 * @param row	the row's index to change
	 * @param vec	the new row vector
	 */
	public void setRow(int row, DenseVector vec) {
		for (int j = 0; j < cols; ++j) {
			data[row][j] = vec.get(j);
		}
	}
	
	/**
	 * Returns the opposite dense matrix.
	 * 
	 * @return the opposite dense matrix
	 */
	public DenseMatrix opposite() {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = -data[i][j];
			}
		}
		return result;
	}
	
	/**
	 * <p>Addition between dense matrices.</p>
	 * Adds the matrix {@code that}.
	 * 
	 * @param that	the matrix to add
	 * @return the resulting matrix
	 */
	public DenseMatrix add(DenseMatrix that) {
		assert(rows == that.rows && cols == that.cols);
		
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] + that.data[i][j];
			}
		}
		return result;
	}
	
	/**
	 * Adds the given {@code value} to all the entries.
	 * 
	 * @param value	 the value to add
	 * @return the resulting matrix
	 */
	public DenseMatrix add(double value) {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] + value;
			}
		}
		return result;
	}
	
	/**
	 * Adds the given {@code value} to the entry in position
	 * ({@code row}, {@code col}).
	 * 
	 * @param row	 the entry's row
	 * @param col	 the entry's column
	 * @param value	 the value to add
	 */
	public void add(int row, int col, double value) {
		data[row][col] += value;
	}
	
	/**
	 * <p>Subtraction between dense matrices.</p>
	 * Subtracts the matrix {@code that}.
	 * 
	 * @param that	the matrix to subtract
	 * @return the resulting matrix
	 */
	public DenseMatrix sub(DenseMatrix that) {
		assert(rows == that.rows && cols == that.cols);
		
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] - that.data[i][j];
			}
		}
		return result;
	}
	
	/**
	 * Subtracts the given {@code value} to all the entries.
	 * 
	 * @param value	 the value to subtract
	 * @return the resulting matrix
	 */
	public DenseMatrix sub(double value) {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] - value;
			}
		}
		return result;
	}
	
	/**
	 * Subtracts the given {@code value} to the entry in position
	 * ({@code row}, {@code col}).
	 * 
	 * @param row	 the entry's row
	 * @param col	 the entry's column
	 * @param value	 the value to subtract
	 */
	public void sub(int row, int col, double value) {
		data[row][col] -= value;
	}
	
	/**
	 * <p>Dot product of dense matrices.</p>
	 * Applies the dot product between this and {@code that} matrix.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public DenseMatrix dot(DenseMatrix that) {
		assert(cols == that.rows);
		
		DenseMatrix result = new DenseMatrix(rows, that.cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < that.cols; ++j) {
				for (int k = 0; k < cols; ++k) {
					result.data[i][j] += data[i][k] * that.data[k][j]; 
				}
			}
		}
		return result;
	}

	/**
	 * <p>Point-wise multiplication of dense matrices.</p>
	 * Multiplies point-wise this matrix by the matrix {@code that}.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public DenseMatrix hadamard(DenseMatrix that) {
		assert(rows == that.rows && cols == that.cols);

		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] * that.data[i][j];
			}
		}
		return result;
	}

	/**
	 * <p>Product between matrix and vector.</p>
	 * Applies the product between this vector and the given {@code matrix}. </br>
	 * The product between an <tt>n x m</tt> matrix
	 * and a <tt>m x 1</tt> matrix (i.e., n-dimensional vector) is an
	 * <tt>n x 1</tt> matrix (i.e., column vector). 
	 * 
	 * <p>See also {@link DenseMatrix#times(DenseVector)}.</p>
	 * 
	 * @param vec	the vector
	 * @return the resulting vector
	 */
	public DenseVector times(DenseVector vec) {
		assert(cols == vec.size());

		DenseVector result = new DenseVector(rows);
		for (int i = 0; i < rows; i++) {
			result.set(i, row(i).dot(vec));
		}
		return result;
	}
	
	/**
	 * Computes the transposition of the matrix.
	 * 
	 * @return the transposition of the matrix
	 */
	public DenseMatrix transpose() {
		DenseMatrix result = new DenseMatrix(cols, rows);
		for (int i = 0; i < result.rows; i++) {
			for (int j = 0; j < result.cols; j++) {
				result.data[i][j] =  data[j][i];
			}
		}
		return result;
	}
	
	/**
	 * Computes the trace of the matrix: sum of all values
	 * in the main diagonal.
	 * 
	 * @return the trace of the matrix
	 */
	public double trace() {
		assert(isSquare());
		
		double result = 0;
		for (int i = 0; i < rows; ++i) {
			result += data[i][i];
		}
		return result;
	}
	
	/**
	 * <p>Computes the inverse of the matrix.</p>
	 * <strong>This implementation is adopted from PREA package.</strong>
	 * 
	 * @return the inverse matrix
	 */
	public DenseMatrix inv() {
		assert(isSquare());

		DenseMatrix result = DenseMatrix.identity(rows);

		if (rows == 1) {
			result.set(0, 0, 1 / data[0][0]);
			return result;
		}

		DenseMatrix b = new DenseMatrix(this);
		for (int i = 0; i < rows; ++i) {
			// find pivot:
			double mag = 0;
			int pivot = -1;

			for (int j = i; j < rows; ++j) {
				double mag2 = Math.abs(b.get(j, i));
				if (mag2 > mag) {
					mag = mag2;
					pivot = j;
				}
			}

			// no pivot (error):
			if (pivot == -1 || mag == 0) {
				return result;
			}

			// move pivot row into position:
			if (pivot != i) {
				double temp;
				for (int j = i; j < rows; ++j) {
					temp = b.get(i, j);
					b.set(i, j, b.get(pivot, j));
					b.set(pivot, j, temp);
				}

				for (int j = 0; j < rows; ++j) {
					temp = result.get(i, j);
					result.set(i, j, result.get(pivot, j));
					result.set(pivot, j, temp);
				}
			}

			// normalize pivot row:
			mag = b.get(i, i);
			for (int j = i; j < rows; ++j) {
				b.set(i, j, b.get(i, j) / mag);
			}

			for (int j = 0; j < rows; ++j) {
				result.set(i, j, result.get(i, j) / mag);
			}

			// eliminate pivot row component from other rows:
			for (int k = 0; k < rows; ++k) {
				if (k == i) continue;

				double mag2 = b.get(k, i);

				for (int j = i; j < rows; ++j) {
					b.set(k, j, b.get(k, j) - mag2 * b.get(i, j));
				}

				for (int j = 0; j < rows; ++j) {
					result.set(k, j, result.get(k, j) - mag2 * result.get(i, j));
				}
			}
		}

		return result;
	}
	
	/**
	 * Computes the {@code n}-norm of the matrix.
	 * 
	 * @param n		the norm's index
	 * @return the {@code n}-norm of the matrix
	 */
	public double norm(int n) {
		assert(n > 0);
		
		double result = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result += Math.abs(Math.pow(data[i][j], n));
			}
		}
		result = Math.pow(result, 1.0 / (double) n);
		return result;
	}
	
	/**
	 * Converts the matrix to a bi-dimensional vector of doubles.
	 * 
	 * @return the bi-dimensional vector representation of the matrix
	 */
	public double[][] toArray() {
		return data.clone();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		
		if (!(obj instanceof DenseMatrix)) {
			return false;
		}
		
		DenseMatrix that = (DenseMatrix) obj;
		
		if (rows != that.rows || cols != that.cols) {
			return false;
		}
		
		boolean equal = true;
		for (int i = 0; i < rows && equal; ++i) {
			for (int j = 0; j < cols && equal; ++j) {
				equal &= Misc.isEqual(data[i][j], that.data[i][j]);
			}
		}
		
		return equal;
	}
	
	@Override
	public int hashCode() {
		int hash = rows;
		hash = 17 * hash + cols; 
		hash = 31 * hash + data.hashCode();
		return hash;
	}
	
	/**
	 * Converts the matrix to a sparse representation.
	 * 
	 * @return the sparse representation of the matrix
	 */
	public SparseMatrix toSparse() {
		ArrayList<Integer> r = new ArrayList<Integer>();
		ArrayList<Integer> c = new ArrayList<Integer>();
		ArrayList<Double> v = new ArrayList<Double>();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (!Misc.isEqual(data[i][j], 0.0)) {
					r.add(i);
					c.add(j);
					v.add(data[i][j]);
				}
			}
		}
		return new SparseMatrix(rows, cols, Ints.toArray(r), Ints.toArray(r), Doubles.toArray(v)); 
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Dimension: " + rows + " x " + cols + "\n");

		for (int i = 0; i < rows; ++i) {
			sb.append("[");
			for (int j = 0; j < cols; ++j) {
				sb.append((float) data[i][j]);
				if (j < cols - 1) {
					sb.append("\t");
				}
			}
			sb.append("]\n");
		}
		return sb.toString();
	}
}
