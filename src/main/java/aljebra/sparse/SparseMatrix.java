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
package aljebra.sparse;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import aljebra.dense.DenseMatrix;
import aljebra.dense.DenseVector;
import aljebra.utils.Misc;
import aljebra.utils.RNG;

/**
 * <p>Sparse implementation of a matrix.</p>
 * <p>The implementation is suitable for <i>immutable</i> sparse matrix, the 
 * usage of this class as a dense matrix can have bad performances with respect to
 * the dense implementation provided by the class {@link DenseMatrix}.</p>
 * 
 * <strong>This implementation is a re-elaboration of the one
 * provided by Guibing Guo in his <a href="http://www.librec.net">librec</a>.</strong>
 * 
 * @author Mirko Polato
 *
 */
public class SparseMatrix implements Serializable, Cloneable {

	private static final long serialVersionUID = -3531663811502120771L;
	
	//
	// STATIC SECTION
	//
	
	/**
	 * Constructs a new zero matrix with dimension 
	 * {@code rows} x {@code cols}.
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 * @return a new sparse null matrix
	 */
	public static SparseMatrix zero(int rows, int cols) {
		return new SparseMatrix(rows, cols, true);
	}
	
	/**
	 * <p>Constructs a new matrix with all entries set to one,
	 * with dimension {@code rows} x {@code cols}.</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse matrix
	 * is not as efficient as using a dense implementation. See {@link DenseMatrix}.</strong>
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 * @return a new one matrix
	 */
	public static SparseMatrix one(int rows, int cols) {
		SparseMatrix one = new SparseMatrix(rows, cols, false);
		int size = rows * cols;
		
		one.rowPtr = new int[rows + 1];
		one.colInd = new int[size];
		one.rowData = new double[size];

		one.colPtr = new int[cols + 1];
		one.rowInd = new int[size];
		one.colData = new double[size];
		
		one.rowPtr[0] = one.colPtr[0] = 0;
		for (int i = 0; i < rows; ++i) {
			one.rowPtr[i + 1] = (i + 1) * cols;
			for (int j = 0; j < cols; ++j) {
				int k = i * cols + j;
				one.colInd[k] = j;
				one.rowData[k] = 1.0;
			}
		}
		
		for (int j = 0; j < cols; ++j) {
			one.colPtr[j + 1] = (j + 1) * rows;
			for (int i = 0; i < rows; ++i) {
				int k = j * rows + i;
				one.rowInd[k] = i;
				one.colData[k] = 1.0;
			}
		}
		return one;
	}
	
	/**
	 * <p>Constructs a new sparse matrix initialized with random
	 * double values in the range [0,1].</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse matrix
	 * is not as efficient as using a dense implementation. See {@link DenseMatrix}.</strong>
	 * 
	 * @param rows	number of rows
	 * @param cols	number of columns
	 * @return a new sparse matrix.
	 */
	public static SparseMatrix rand(int rows, int cols) {
		int[] r = new int[rows * cols];
		int[] c = new int[rows * cols];
		double[] v = new double[rows * cols];
		for (int i = 0; i < rows * cols; ++i) {
			for (int j = 0; j < cols; ++j) {
				int k = i * cols + j;
				double val = RNG.uniformDbl();
				r[k] = i;
				c[k] = j;
				v[k] = val;
			}
		}
		clean(r, c, v);
		return new SparseMatrix(r, c, v);
	}
	
	/**
	 * Constructs an identity matrix with dimension 
	 * {@code n} x {@code n}.
	 * 
	 * @param n	 number of rows and column
	 * @return a new identity matrix
	 */
	public static SparseMatrix identity(int n) {
		SparseMatrix eye = new SparseMatrix(n, n, false);
		
		eye.rowPtr = new int[n + 1];
		eye.colInd = new int[n];
		eye.rowData = new double[n];

		eye.colPtr = new int[n + 1];
		eye.rowInd = new int[n];
		eye.colData = new double[n];
		
		eye.rowPtr[0] = eye.colPtr[0] = 0;
		for (int i = 0; i < n; ++i) {
			eye.rowPtr[i + 1] = eye.colPtr[i + 1] = i;
			eye.rowInd[i] = eye.colInd[i] = i;
			eye.rowData[i] = eye.colData[i] = 1.0;
		}
		return eye;
	}
	
	// Clean up the vectors: remove the entry corresponding to zero values.
	private static void clean(int[] r, int[] c, double[] v) {
		assert (r.length == c.length && 
				r.length == v.length);
				
		ArrayList<Integer> skip = new ArrayList<Integer>();
		for (int i = 0; i < v.length; ++i) {
			if (v[i] == 0.0) {
				skip.add(i);
			}
		}
		
		if (skip.size() > 0) {
			int size = r.length - skip.size();
			
			int[] newr = new int[size];
			int[] newc = new int[size];
			double[] newv = new double[size];
			
			for (int i = 0, j = 0; j < skip.size(); ++j) {
				for (int k = i + j; k < skip.get(j); ++k, ++i) {
					newr[i] = r[k];
					newc[i] = c[k];
					newv[i] = v[k];
				}
			}
			
			r = newr;
			c = newc;
			v = newv;
		}
	}
	
	//
	// END STATIC SECTION
	//
	
	
	// Number of rows
	protected final int rows;
	
	// Number of columns
	protected final int cols;
	
	// Compressed Row Storage (CRS)
	protected double[] rowData;
	protected int[] rowPtr, colInd;

	// Compressed Column Storage (CCS)
	protected double[] colData;
	protected int[] colPtr, rowInd;
	
	// Private constructor useful in other methods
	private SparseMatrix(int rows, int cols, boolean init) {
		this.rows = rows;
		this.cols = cols;
		
		if (init) {
			rowPtr = new int[rows + 1];
			colInd = new int[0];
			rowData = new double[0];

			colPtr = new int[cols + 1];
			rowInd = new int[0];
			colData = new double[0];
		}
	}
	
	/**
	 * Constructs a new sparse vector using the given coordinates and
	 * non-zero values. The dimensions of the matrix are inferred
	 * by the given parameters.
	 * 
	 * @param rows		the list of non-zero row indexes
	 * @param cols		the list of non-zero column indexes
	 * @param values	the list of non-zero values
	 */
	public SparseMatrix(int[] rows, int[] cols, double[] values) {
		this(Misc.max(rows) + 1, Misc.max(cols) + 1, rows, cols, values);
	}
	
	/**
	 * Constructs a new sparse vector using the given coordinates and
	 * non-zero values
	 * 
	 * @param nrows		number of rows
	 * @param ncols		number of columns
	 * @param rows		the list of non-zero row indexes
	 * @param cols		the list of non-zero column indexes
	 * @param values	the list of non-zero values
	 */
	public SparseMatrix(int nrows, int ncols, int[] rows, int cols[], double[] values) {
		
		clean(rows, cols, values);
		
		if (nrows < 0) {
			nrows = Misc.max(rows) + 1;
			ncols = Misc.max(cols) + 1;
		} else {
			assert (nrows >= Misc.max(rows) + 1 &&
					ncols >= Misc.max(cols) + 1);
		}
		
		this.rows = nrows;
		this.cols = ncols;
		
		init(rows, cols, values);
	}
	
	// Initializes the structures used to construct the sparse matrix
	private void init(int[] r, int[] c, double[] v) {
		int nnz = v.length;

		rowPtr = new int[rows + 1];
		colInd = new int[nnz];
		rowData = new double[nnz];

		colPtr = new int[cols + 1];
		rowInd = new int[nnz];
		colData = new double[nnz];
		
		rowPtr[0] = colPtr[0] = 0;
		
		Misc.sort(r, c, v);
		
		int j = 0;
		for (int i = 0; i < rows - 1; ++i) {
			while (r[j] == i) {
				++j;
			}
			rowPtr[i + 1] = j;
		}
		rowPtr[rows] = nnz;
		
		for (int i = 1; i <= rows; ++i) {
			Misc.sort(rowPtr[i - 1], rowPtr[i] - 1, c, null, v);
		}
		
		for (int i = 0; i < nnz; ++i) {
			colInd[i] = c[i];
			rowData[i] = v[i];
		}

		Misc.sort(c, r, v);
		
		j = 0;
		for (int i = 0; i < cols - 1; ++i) {
			while (c[j] == i) {
				++j;
			}
			colPtr[i + 1] = j;
		}
		colPtr[cols] = nnz;
		
		for (int i = 0; i < nnz; ++i) {
			rowInd[i] = r[i];
			colData[i] = v[i];
		}
		
		for (int i = 1; i <= cols; ++i) {
			Misc.sort(colPtr[i - 1], colPtr[i] - 1, rowInd, null, colData);
		}
	}
	
	/**
	 * Copy constructors.
	 * 
	 * @param mat	the matrix
	 */
	public SparseMatrix(SparseMatrix mat) {
		this.rows = mat.rows;
		this.cols = mat.cols;

		copyCRS(mat.rowData, mat.rowPtr, mat.colInd);
		copyCCS(mat.colData, mat.colPtr, mat.rowInd);
	}

	// Copy the structures of the "Compressed Row Storage" format
	private void copyCRS(double[] data, int[] ptr, int[] idx) {
		rowData = data.clone();
		rowPtr = ptr.clone();
		colInd = idx.clone();
	}

	// Copy the structures of the "Compressed Column Storage" format
	private void copyCCS(double[] data, int[] ptr, int[] idx) {
		colData = data.clone();
		colPtr = ptr.clone();
		rowInd = idx.clone();
	}

	@Override
	public SparseMatrix clone() {
		return new SparseMatrix(this);
	}
	
	/**
	 * Returns the count of non-zero entries in the matrix.
	 * 
	 * @return the count of non-zero entries in the matrix
	 */
	public int nnzCount() {
		return rowData.length;
	}
	
	/**
	 * Returns the number of rows.
	 * 
	 * @return the number of rows
	 */
	public int rows() {
		return rows;
	}
	
	/**
	 * Returns the number of columns.
	 * 
	 * @return the number of columns
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
	 * Returns the list of row indexes in the "Coordinate" format.
	 * 
	 * @return the list of row indexes in the "Coordinate" format
	 */
	public int[] getCOORows() {
		int[] result = new int[nnzCount()];
		int j = 0;
		for (int i = 0; i < rowPtr.length; ++i) {
			for (int k = rowPtr[i + 1]; k < rowPtr[i]; ++k, ++j) {
				result[j] = i;
			}
		}
		return result;
	}
	
	/**
	 * Returns the list of column indexes in the "Coordinate" format.
	 * 
	 * @return the list of column indexes in the "Coordinate" format
	 */
	public int[] getCOOCols() {
		return colInd.clone();
	}
	
	/**
	 * Returns the list of values in the "Coordinate" format.
	 * 
	 * @return the list of values in the "Coordinate" format
	 */
	public double[] getCOOValues() {
		return rowData.clone();
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
		int index = Arrays.binarySearch(colInd, rowPtr[row], rowPtr[row + 1], col);
		return (index >= 0) ? rowData[index] : 0.0;
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
		
		int i = Arrays.binarySearch(colInd, rowPtr[row], rowPtr[row + 1], col);
		int j = Arrays.binarySearch(rowInd, colPtr[col], colPtr[col + 1], row);
		
		int nnz = rowData.length;
		boolean reduce = false;

		if (i >= 0) {
			if (!Misc.isEqual(value, 0.0)) {
				rowData[i] = value;
				colData[j] = value;
				return;
			} else {
				--nnz;
				reduce = true;
			}
		} else {
			if (Misc.isEqual(value, 0.0)) {
				return;
			}
			
			++nnz;
			i = -(i + 1);
			j = -(j + 1);
		}
		
		int[] newRowInd = new int[nnz];
		int[] newColInd = new int[nnz];
		double[] newRowData = new double[nnz];
		double[] newColData = new double[nnz];

		System.arraycopy(rowInd, 0, newRowInd, 0, j);
		System.arraycopy(colInd, 0, newColInd, 0, i);
		System.arraycopy(rowData, 0, newRowData, 0, i);
		System.arraycopy(colData, 0, newColData, 0, j);

		int delta = 1;
		if (reduce) {
			++i;
			++j;
			delta = -1;
		} else {
			newRowInd[j] = row;
			newColInd[i] = col;
			newRowData[i] = newColData[j] = value;
		}
		
		System.arraycopy(rowInd, j, newRowInd, j + delta, nnz - j - delta);
		System.arraycopy(colInd, i, newColInd, i + delta, nnz - i - delta);
		System.arraycopy(rowData, i, newRowData, i + delta, nnz - i - delta);
		System.arraycopy(colData, j, newColData, j + delta, nnz - j - delta);

		rowInd = newRowInd;
		colInd = newColInd;
		rowData = newRowData;
		colData = newColData;
		
		for (int k = row + 1; k <= rows; ++k) {
			rowPtr[k] += delta;
		}
		
		for (int k = col + 1; k <= cols; ++k) {
			colPtr[k] += delta;
		}
	}
	
	/**
	 * Returns the {@code i}-th row as a sparse vector.
	 * 
	 * @param i	 the row's index
	 * @return the {@code i}-th row as a sparse vector
	 */
	public SparseVector row(int i) {
		assert(i >= 0 && i < rows);
		
		SparseVector sv = new SparseVector(cols);
		for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
			sv.set(colInd[j], rowData[j]);
		}
		return sv;
	}
	
	/**
	 * Returns the {@code j}-th column as a sparse vector.
	 * 
	 * @param j  the column's index
	 * @return the {@code j}-th column as a sparse vector
	 */
	public SparseVector col(int j) {
		assert(j >= 0 && j < cols);
		
		SparseVector sv = new SparseVector(rows);
		for (int i = colPtr[j]; i < rowPtr[j + 1]; ++i) {
			sv.set(rowInd[i], colData[i]);
		}
		return sv;
	}
	
	/**
	 * Returns the main diagonal of the matrix as a 
	 * sparse vector.
	 * 
	 * @return the main diagonal of the matrix as a sparse vector
	 */
	public SparseVector diag() {
		assert(isSquare());
		
		SparseVector result = new SparseVector(rows);
		for (int i = 0; i < rows; ++i) {
			for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
				if (colInd[j] == i) {
					result.set(i, rowData[j]);
					break;
				}
			}
		}
		return result;
	}
	
	/**
	 * Returns the opposite sparse matrix.
	 * 
	 * @return the opposite sparse matrix
	 */
	public SparseMatrix opposite() {
		SparseMatrix result = new SparseMatrix(rows, cols, true);
		
		double[] rowDataNeg = new double[rowData.length];
		double[] colDataNeg = new double[colData.length];
		
		for (int i = 0; i < rowData.length; ++i) {
			rowDataNeg[i] = -rowData[i];
		}
		for (int i = 0; i < colData.length; ++i) {
			colDataNeg[i] = -colData[i];
		}
		
		result.copyCRS(rowDataNeg, rowPtr, colInd);
		result.copyCCS(colDataNeg, colPtr, rowInd);
		
		return result;
	}
	
	/**
	 * <p>Addition between sparse matrices.</p>
	 * Adds the matrix {@code that}.
	 * 
	 * @param that	the matrix to add
	 * @return the resulting matrix
	 */
	public SparseMatrix add(SparseMatrix that) {
		assert(rows == that.rows && cols == that.cols);
		
		SparseMatrix result = new SparseMatrix(rows, cols, true);
		for (int i = 0; i < rows; ++i) {
			int j = rowPtr[i];
			int k = that.rowPtr[i];
			while (j < rowPtr[i + 1] || k < that.rowPtr[i + 1]) {
				int cj = colInd[j];
				int ck = that.colInd[k];
				if (cj < ck) {
					result.set(i, cj, rowData[j++]);
				} else if (cj == ck) {
					result.set(i, cj, rowData[j++] + that.rowData[k++]);
				} else {
					result.set(i, ck, that.rowData[k++]);
				}
			}
		}
		return result;
	}
	
	/**
	 * <p>Subtraction between sparse matrices.</p>
	 * Subtracts the matrix {@code that}.
	 * 
	 * @param that	the matrix to subtract
	 * @return the resulting matrix
	 */
	public SparseMatrix sub(SparseMatrix that) {
		assert(rows == that.rows && cols == that.cols);
		
		SparseMatrix result = new SparseMatrix(rows, cols, true);
		for (int i = 0; i < rows; ++i) {
			int j = rowPtr[i];
			int k = that.rowPtr[i];
			while (j < rowPtr[i + 1] || k < that.rowPtr[i + 1]) {
				int cj = colInd[j];
				int ck = that.colInd[k];
				if (cj < ck) {
					result.set(i, cj, -rowData[j++]);
				} else if (cj == ck) {
					result.set(i, cj, rowData[j++] - that.rowData[k++]);
				} else {
					result.set(i, ck, -that.rowData[k++]);
				}
			}
		}
		return result;
	}
	
	/**
	 * <p>Point-wise multiplication of sparse matrices.</p>
	 * Multiplies point-wise this matrix by the matrix {@code that}.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public SparseMatrix hadamard(SparseMatrix that) {
		assert(rows == that.rows && cols == that.cols);
		
		SparseMatrix result = new SparseMatrix(rows, cols, true);
		for (int i = 0; i < rows; ++i) {
			int j = rowPtr[i];
			int k = that.rowPtr[i];
			while (j < rowPtr[i + 1] || k < that.rowPtr[i + 1]) {
				int cj = colInd[j];
				int ck = that.colInd[k];
				if (cj < ck) {
					++j;
				} else if (cj == ck) {
					result.set(i, cj, rowData[j++] * that.rowData[k++]);
				} else {
					++k;
				}
			}
		}
		return result;
	}
	
	/**
	 * <p>Dot product of sparse matrices.</p>
	 * Applies the dot product between this and {@code that} matrix.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public SparseMatrix dot(SparseMatrix that) {
		assert(rows == that.cols && cols == that.rows);
		
		SparseMatrix result = new SparseMatrix(rows, that.cols, true);
		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < cols; ++c) {
				double res = 0;
				
				int j = rowPtr[r];
				int k = that.colPtr[c];
				while (j < rowPtr[r + 1] || k < that.colPtr[c + 1]) {
					int cj = colInd[j];
					int rk = that.rowInd[k];
					if (cj < rk) {
						++j;
					} else if (cj == rk) {
						res += rowData[j++] * that.colData[k++];
					} else {
						++k;
					}
				}
				
				result.set(r, c, res);
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
	public SparseVector times(SparseVector vec) {
		assert(cols == vec.size());
		
		SparseVector result = new SparseVector(rows);
		for (int i = 0; i < rows; ++i) {
			double res = 0;
			for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
				res += rowData[j] * vec.get(colInd[j]);
			}
			result.set(i, res);
		}
		return result;
	}
	
	/**
	 * Computes the transposition of the matrix.
	 * 
	 * @return the transposition of the matrix
	 */
	public SparseMatrix transpose() {
		SparseMatrix tr = new SparseMatrix(cols, rows, false);

		tr.copyCCS(this.rowData, this.rowPtr, this.colInd);
		tr.copyCRS(this.colData, this.colPtr, this.rowInd);

		return tr;
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
			for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
				if (colInd[j] == i) {
					result += rowData[j];
					break;
				}
			}
		}
		return result;
	}
	
	/**
	 * Computes the sum of all entries.
	 * 
	 * @return the sum of all entries
	 */
	public double sum() {
		double sum = 0;
		for (double d : rowData) {
			sum += d;
		}
		return sum;
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
			for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
				result += Math.abs(Math.pow(rowData[j], n));
			}
		}
		result = Math.pow(result, 1.0 / (double) n);
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		
		if (!(obj instanceof SparseMatrix)) {
			return false;
		}
		
		SparseMatrix that = (SparseMatrix) obj;
		
		if (rows != that.rows || cols != that.cols) {
			return false;
		}
		
		boolean equal = true;
		for (int i = 0; i < rows && equal; ++i) {
			for (int j = rowPtr[i]; j < rowPtr[i + 1] && equal; ++j) {
				double x = rowData[j];
				equal &= Misc.isEqual(x, that.get(i, colInd[j]));
			}
		}
		return equal;
	}
	
	@Override
	public int hashCode() {
		int hash = rows;
		hash = 7 * hash + cols; 
		hash = 13 * hash + rowPtr.hashCode();
		hash = 17 * hash + colInd.hashCode();
		hash = 31 * hash + rowData.hashCode();
		return hash;
	}
	
	/**
	 * Converts the sparse matrix to an array of doubles.
	 * 
	 * @return the corresponding array of doubles
	 */
	public double[][] toArray() {
		double[][] result = new double[rows][cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
				result[i][colInd[j]] = rowData[j];
			}
		}
		return result;
	}
	
	/**
	 * Converts the sparse matrix to a dense representation.
	 * 
	 * @return the dense representation of the matrix
	 */
	public DenseMatrix toDense() {
		DenseMatrix result = new DenseMatrix(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
				result.set(i, colInd[j], rowData[j]);
			}
		}
		return result;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for (int i = 0; i < rows; i++) {
			if (i == 0) {
				sb.append("[");
			} else {
				sb.append(" [");
			}
			for (int j = 0; j < cols; j++) {
				sb.append(get(i, j));
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
}
