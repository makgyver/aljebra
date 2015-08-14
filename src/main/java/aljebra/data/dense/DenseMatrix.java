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
package aljebra.data.dense;

import java.util.ArrayList;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import aljebra.data.IFunction;
import aljebra.data.IMatrix;
import aljebra.data.IVector;
import aljebra.data.sparse.SparseMatrix;
import aljebra.lapack.SVD;
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
public class DenseMatrix implements IMatrix {

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
		assert(M.cols == N.rows);
		assert(mrow >= 0 && mrow < M.rows && ncol >= 0 && ncol < N.cols);

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
	
	@Override
	public DenseMatrix copy() {
		return new DenseMatrix(this);
	}

	@Override
	public double get(int row, int col) {
		assert(row >= 0 && row < rows && col >= 0 && col < cols);
		return data[row][col];
	}
	
	@Override
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
	
	@Override
	public int rows() {
		return rows;
	}
	
	@Override
	public int cols() {
		return cols;
	}
	
	@Override
	public boolean isSquare() {
		return rows == cols;
	}
	
	@Override
	public DenseVector row(int i) {
		assert(i >= 0 && i < rows);
		return new DenseVector(data[i]);
	}
	
	@Override
	public DenseVector col(int j) {
		assert(j >= 0 && j < cols);
		
		DenseVector result = new DenseVector(rows);
		for (int i = 0; i < rows; ++i) {
			result.set(i, data[i][j]);
		}
		return result;
	}
	
	@Override
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
	public void setRow(int row, IVector vec) {
		assert(row >= 0 && row < rows);
		
		for (int j = 0; j < cols; ++j) {
			data[row][j] = vec.get(j);
		}
	}
	
	@Override
	public DenseMatrix opposite() {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = -data[i][j];
			}
		}
		return result;
	}
	
	@Override
	public DenseMatrix add(IMatrix that) {
		assert(rows == that.rows() && cols == that.cols());
		
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] + that.get(i, j);
			}
		}
		return result;
	}
	
	@Override
	public DenseMatrix add(double value) {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] + value;
			}
		}
		return result;
	}
	
	@Override
	public void add(int row, int col, double value) {
		assert(row >= 0 && row < rows && col >= 0 && col < cols);
		data[row][col] += value;
	}
	
	@Override
	public DenseMatrix sub(IMatrix that) {
		assert(rows == that.rows() && cols == that.cols());
		
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] - that.get(i, j);
			}
		}
		return result;
	}
	
	@Override
	public DenseMatrix sub(double value) {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] - value;
			}
		}
		return result;
	}
	
	@Override
	public void sub(int row, int col, double value) {
		assert(row >= 0 && row < rows && col >= 0 && col < cols);
		data[row][col] -= value;
	}
	
	@Override
	public DenseMatrix dot(IMatrix that) {
		assert(cols == that.rows());
		
		DenseMatrix result = new DenseMatrix(rows, that.cols());
		for (int i = 0; i < rows; ++i) {
			for (int k = 0; k < cols; ++k) {
				for (int j = 0; j < that.cols(); ++j) {
					result.data[i][j] += data[i][k] * that.get(k ,j); 
				}
			}
		}
		return result;
	}

	@Override
	public DenseMatrix hadamard(IMatrix that) {
		assert(rows == that.rows() && cols == that.cols());

		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] * that.get(i, j);
			}
		}
		return result;
	}

	@Override
	public DenseVector times(IVector vec) {
		assert(cols == vec.size());

		DenseVector result = new DenseVector(rows);
		for (int i = 0; i < rows; ++i) {
			result.set(i, row(i).dot(vec));
		}
		return result;
	}
	
	@Override
	public DenseMatrix scale(double factor) {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = factor * data[i][j];
			}
		}
		return result;
	}
	
	@Override
	public DenseMatrix div(IMatrix that) {
		assert(rows == that.rows() && cols == that.cols());

		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (that.get(i, j) != 0) {
					result.data[i][j] = data[i][j] / that.get(i, j);
				}
			}
		}
		return result;
	}

	@Override
	public DenseMatrix pow(double power) {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = Math.pow(data[i][j], power);
			}
		}
		return result;
	}
	
	@Override
	public DenseMatrix transpose() {
		DenseMatrix result = new DenseMatrix(cols, rows);
		for (int i = 0; i < result.rows; i++) {
			for (int j = 0; j < result.cols; j++) {
				result.data[i][j] =  data[j][i];
			}
		}
		return result;
	}
	
	@Override
	public double trace() {
		assert(isSquare());
		
		double result = 0;
		for (int i = 0; i < rows; ++i) {
			result += data[i][i];
		}
		return result;
	}
	
	@Override
	public double sum() {
		double sum = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				sum += data[i][j];
			}
		}
		return sum;
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
	
	@Override
	public double mean() {
		return sum() / (rows * cols);
	}
	
	@Override
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
	 * <p>Performs the singular value decomposition of the matrix.</p>
	 * For more information about SVD see the 
	 * <a href="https://en.wikipedia.org/wiki/Singular_value_decomposition">wikipedia page</a>.
	 * 
	 * @return the SVD decomposition object
	 */
	public SVD svd() {
		SVD svd = new SVD(this);
		return svd;
	}
	
	@Override
	public double[][] toArray() {
		double[][] array = new double[rows][cols];
		for (int i = 0; i < rows; ++i) {
			array[i] = data[i].clone();
		}
		return array;
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
		sb.append("[");
		for (int i = 0; i < rows; ++i) {
			if (i == 0) {
				sb.append("[");
			} else {
				sb.append(" [");
			}
			for (int j = 0; j < cols; ++j) {
				sb.append((float) data[i][j]);
				if (j < cols - 1) {
					sb.append("\t");
				}
			}
			sb.append("]");
			if (i == rows - 1) {
				sb.append("]");
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	@Override
	public void apply(IFunction fun) {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				data[i][j] = fun.apply(data[i][j]);
			}
		}
	}
}
