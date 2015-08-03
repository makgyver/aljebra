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

import aljebra.dense.DenseVector;
import aljebra.utils.Misc;
import aljebra.utils.RNG;

/**
 * <p>Sparse implementation of a vector.</p>
 * <p>The implementation is suitable for <i>immutable</i> sparse vector, the 
 * usage of this class as a dense vector can have bad performances with respect to
 * the dense implementation provided by the class {@link DenseVector}.</p>
 * 
 * <strong>This implementation is a re-elaboration of the one
 * provided by Guibing Guo in his <a href="http://www.librec.net">librec</a>.</strong>
 * 
 * @author Mirko Polato
 *
 */
public class SparseVector implements Serializable, Cloneable {

	private static final long serialVersionUID = 849869514972550140L;
	
	//
	// STATIC SECTION
	//
	
	/**
	 * Constructs a new sparse vector with all entries set to 0.
	 * 
	 * @param size	the size of the vector
	 * @return a new sparse null vector
	 */
	public static SparseVector zero(int size) {
		return new SparseVector(size);
	}
	
	/**
	 * <p>Constructs a new sparse vector with all entries set to 1.</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse vector
	 * is not as efficient as using a dense implementation. See {@link DenseVector}.</strong>
	 * 
	 * @param size 	the size of the vector
	 * @return a new sparse one vector
	 */
	public static SparseVector one(int size) {
		SparseVector result = new SparseVector(size);
		result.setAll(1.0);
		return result;
	}
	
	/**
	 * <p>Constructs a new sparse vector initialized with random
	 * double values in the range [0,1].</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse vector
	 * is not as efficient as using a dense implementation. See {@link DenseVector}.</strong>
	 * 
	 * @param size	the size of the vector
	 * @return a new sparse vector.
	 */
	public static SparseVector rand(int size) {
		int[] ids = new int[size];
		double[] vals = new double[size];
		for (int i = 0; i < size; ++i) {
			ids[i] = i;
			vals[i] = RNG.uniformDbl();
		}
		return new SparseVector(ids, vals);
	}
	
	//
	// END STATIC SECTION
	//
	
	
	// The actual vector's data (only non-zero values)
	protected double[] data;
	
	// The indexes of the corresponding values in 'data'
	protected int[] ids;
	
	// Counts the number of non-zero values in the vector
	protected int count;
	
	// The vector's size
	protected final int size;
	
	/**
	 * Creates a new null sparse vector.
	 * 
	 * @param size	the size of the vector
	 */
	public SparseVector(int size) {
		this.size = size;
		this.ids = new int[0];
		this.data = new double[0];
	}
	
	/**
	 * Constructs a new sparse vector with the given non-zero values.
	 * 
	 * @param size		the size of the vector
	 * @param indexes	the list of non-zero indexes
	 * @param values	the list of non-zero values
	 */
	public SparseVector(int size, int[] indexes, double[] values) {
		clean(indexes, values);
		
		assert(indexes.length == values.length && size >= Misc.max(indexes) + 1);
		
		this.size = size;
		this.ids = indexes;
		this.data = values;
		this.count = indexes.length;
	}
	
	/**
	 * Constructs a new sparse vector with the given non-zero values.
	 * The size of the vector is inferred from the parameters.
	 * 
	 * @param indexes	the list of non-zero indexes
	 * @param values	the list of non-zero values
	 */
	public SparseVector(int[] indexes, double[] values) {
		this(Misc.max(indexes) + 1, indexes, values);
	}
	
	/**
	 * Copy constructor.
	 * 
	 * @param vec	the sparse vector to copy
	 */
	public SparseVector(SparseVector vec) {
		this.size = vec.size;
		this.data = new double[size];
		this.ids = new int[size];
		for (int i = 0; i < size; ++i) {
			this.data[i] = vec.data[i];
			this.ids[i] = vec.ids[i];
		}
		this.count = vec.count;
	}
	
	// Clean up the vectors: remove the entry corresponding to zero values.
	private void clean(int[] ids, double[] v) {
		assert (ids.length == v.length);
				
		ArrayList<Integer> skip = new ArrayList<Integer>();
		for (int i = 0; i < v.length; ++i) {
			if (!Misc.isEqual(v[i], 0.0)) {
				skip.add(i);
			}
		}
		
		if (skip.size() > 0) {
			int size = ids.length - skip.size();
			
			int[] newi = new int[size];
			double[] newv = new double[size];
			
			for (int i = 0, j = 0; j < skip.size(); ++j) {
				for (int k = i + j; k < skip.get(j); ++k, ++i) {
					newi[i] = ids[k];
					newv[i] = v[k];
				}
			}
			
			ids = newi;
			v = newv;
		}
	}
	
	/**
	 * Returns the value in position {@code index}.
	 * 
	 * @param index	 the entry's position
	 * @return the value in position {@code index}
	 */
	public double get(int index) {
		assert(index >= 0 && index < size);
		
		int i = Arrays.binarySearch(ids, 0, count, index);
		return i >= 0 ? data[i] : 0;
	}
	
	/**
	 * Sets the entry in position {@code index} to the given {@code value}.
	 * 
	 * @param index	 the entry's position
	 * @param value	 the new value (expected non-zero)
	 */
	public void set(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = getIndex(index);
			data[i] = value;
		}
	}
	
	/**
	 * <p>Sets all entries to the given {@code value}.</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse vector
	 * is not as efficient as using a dense implementation. See {@link DenseVector}.</strong>
	 * 
	 * @param value	 the new value for the entries
	 */
	public void setAll(double value) {
		if (!Misc.isEqual(value, 0.0)) {
			data = new double[size];
			ids = new int[size];
			for (int i = 0; i < size; ++i) {
				data[i] = value;
				ids[i] = i;
			}
			count = size;
		} else {
			data = new double[0];
			ids = new int[0];
			count = 0;
		}
	}
	
	/**
	 * Returns the size of the vector.
	 * 
	 * @return the size of the vector
	 */
	public int size() {
		return size;
	}
	
	/**
	 * Returns the number of non-zero values in the vector.
	 *  
	 * @return the number of non-zero values in the vector
	 */
	public int nnzCount() {
		return count;
	}
	
	// Tries to find the index. If it is not found, a reallocation is done, 
	// and a new index is returned.
	private int getIndex(int idx) {
 		int i = Arrays.binarySearch(ids, 0, count, idx);

		// Found
		if (i >= 0) {
			return i;
		}

		int[] newIndex = ids;
		double[] newData = data;

		// get insert position
		i = -(i + 1);

		// Check available memory
		if (++count > data.length) {

			// If zero-length, use new length of 1, else double the bandwidth
			int newLength = data.length != 0 ? data.length << 1 : 1;

			// Copy existing data into new arrays
			newIndex = new int[newLength];
			newData = new double[newLength];
			System.arraycopy(ids, 0, newIndex, 0, i);
			System.arraycopy(data, 0, newData, 0, i);
		}

		// All ok, make room for insertion
		System.arraycopy(ids, i, newIndex, i + 1, count - i - 1);
		System.arraycopy(data, i, newData, i + 1, count - i - 1);

		// Put in new structure
		newIndex[i] = idx;
		newData[i] = 0.;

		// Update pointers
		ids = newIndex;
		data = newData;

		// Return insertion index
		return i;
	}
	
	/**
	 * Returns the opposite sparse vector.
	 * 
	 * @return the opposite sparse vector
	 */
	public SparseVector opposite() {
		double[] newData = new double[count];
		for (int i = 0; i < count; ++i) {
			newData[i] = -data[i];
		}
		return new SparseVector(size, ids, newData);
	}
	
	/**
	 * <p>Addition between sparse vectors.</p>
	 * Adds the vector {@code that}.
	 * 
	 * @param that	the vector to add
	 * @return the resulting vector
	 */
	public SparseVector add(SparseVector that) {
		assert(size == that.size);
		
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) + that.get(i));
		}
		return result;
	}
	
	/**
	 * <p>Addition between a sparse vector and a dense one.</p>
	 * Adds the vector {@code that}.
	 * 
	 * @param that	the vector to add
	 * @return the resulting vector
	 */
	public SparseVector add(DenseVector that) {
		assert(size == that.size());
		
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) + that.get(i));
		}
		return result;
	}
	
	
	/**
	 * <p>Adds the given {@code value} to all the vector's entries.</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse vector
	 * is not as efficient as using a dense implementation. See {@link DenseVector}.</strong>
	 * 
	 * @param value	 the value to add
	 * @return the resulting vector
	 */
	public SparseVector add(double value) {
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) + value);
		}
		return result;
	}
	
	/**
	 * Adds the given {@code value} to the entry in position
	 * {@code index}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the value to add
	 */
	public void add(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = getIndex(index);
			data[i] += value;
		}
	}
	
	/**
	 * <p>Subtraction between sparse vectors.</p>
	 * Subtracts the vector {@code that}.
	 * 
	 * @param that	the vector to subtract
	 * @return the resulting vector
	 */
	public SparseVector sub(SparseVector that) {
		assert(size == that.size);
		
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) - that.get(i));
		}
		return result;
	}
	
	/**
	 * <p>Subtraction between a sparse vector and a dense one.</p>
	 * Subtracts the vector {@code that}.
	 * 
	 * @param that	the vector to subtract
	 * @return the resulting vector
	 */
	public SparseVector sub(DenseVector that) {
		assert(size == that.size());
		
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) - that.get(i));
		}
		return result;
	}
	
	/**
	 * <p>Subtract the given {@code value} to all the vector's entries.</p>
	 * 
	 * <strong>NOTE: using a sparse implementation with a non-sparse vector
	 * is not as efficient as using a dense implementation. See {@link DenseVector}.</strong>
	 * 
	 * @param value	 the value to subtract
	 * @return the resulting vector
	 */
	public SparseVector sub(double value) {
		SparseVector result = new SparseVector(size);
		for (int i : ids) {
			result.set(i, data[i] - value);
		}
		return result;
	}
	
	/**
	 * Subtracts the given {@code value} to the entry in position
	 * {@code index}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the value to subtract
	 */
	public void sub(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = getIndex(index);
			data[i] -= value;
		}
	}
	
	/**
	 * Scales the vector by the given {@code factor}.
	 * 
	 * @param factor	the scale factor
	 * @return the scaled vector
	 */
	public SparseVector scale(double factor) {
		if (factor == 0.0) {
			return new SparseVector(size);
		}
		
		SparseVector result = new SparseVector(size);
		for (int i : ids) {
			result.set(i, data[i] * factor);
		}
		return result;
	}
	
	/**
	 * Multiplies the given {@code value} to the entry in position
	 * {@code index}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the multiplicative factor
	 */
	public void times(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = Arrays.binarySearch(ids, index);
			if (i >= 0) {
				data[i] *= value;
			}
		} else {
			data = new double[0];
			ids = new int[0];
			count = 0;
		}
	}
	
	/**
	 * <p>Point-wise multiplication of sparse vectors.</p>
	 * Multiplies point-wise this vector by the vector {@code that}.
	 * 
	 * @param that	the vector to multiply point-wise
	 * @return the resulting vector
	 */
	public SparseVector times(SparseVector that) {
		assert(size == that.size);
		
		SparseVector result = new SparseVector(size);
		for (int i : ids) {
			result.set(i, data[i] * that.data[i]);
		}
		return result;
	}
	
	/**
	 * <p>Product between sparse vector and matrix.</p>
	 * Applies the product between this vector and the given {@code matrix}. </br>
	 * The product between <tt>1 x n</tt> matrix (i.e., n-dimensional vector), 
	 * and an <tt>n x m</tt> matrix is a <tt>1 x m</tt> matrix (i.e., row vector).
	 * 
	 * <p>See also {@link SparseMatrix#times(SparseVector)}.</p>
	 * 
	 * @param matrix	the matrix
	 * @return the resulting vector
	 */
	public SparseVector times(SparseMatrix matrix) {
		assert(size == matrix.rows());
		
		SparseVector result = new SparseVector(matrix.cols());
		for (int i = 0; i < matrix.cols; ++i) {
			double s = 0;
			for (int j : ids) {
				s += data[j] * matrix.get(j, i);
			}
			
			if (s != 0) {
				result.set(i, s);
			}
		}
		return result;
	}
	
	/**
	 * <p>Dot product of sparse vectors.</p>
	 * Applies the dot product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting vector
	 */
	public double dot(SparseVector that) {
		assert(size == that.size);
		
		double result = 0.0;
		for (int i : ids) {
			result += data[i] * that.get(i);
		}
		return result;
	}
	
	/**
	 * <p>Dot product of a sparse vector with a dense one.</p>
	 * Applies the dot product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting vector
	 */
	public double dot(DenseVector that) {
		assert(size == that.size());
		
		double result = 0.0;
		for (int i : ids) {
			result += data[i] * that.get(i);
		}
		return result;
	}
	
	/**
	 * <p>Outer product of sparse vectors.</p>
	 * Applies the outer product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting sparse matrix
	 */
	public SparseMatrix outer(SparseVector that) {
		int[] rows = new int[size * that.size];
		int[] cols = new int[size * that.size];
		double[] vals = new double[size * that.size];
		
		int k = 0;
		for (int i : ids) {
			for (int j : that.ids) {
				rows[k] = i;
				cols[k] = j;
				vals[k] = data[i] * that.data[j];
				++k;
			}
		}
		return new SparseMatrix(size, that.size, rows, cols, vals);
	}
	
	/**
	 * Computes the sum of the vector's entries.
	 * 
	 * @return the sum of the vector's entries
	 */
	public double sum() {
		double result = 0;
		for (double v : data) {
			result += v;
		}
		return result;
	}
	
	/**
	 * Computes the mean of the vector's entries.
	 * 
	 * @return the mean of the vector's entries
	 */
	public double mean() {
		return sum() / size;
	}
	
	/**
	 * Returns the list of non-zero indexes.
	 * 
	 * @return the list of non-zero indexes
	 */
	public int[] indexes() {
		return ids;
	}
	
	/**
	 * Returns the list of non-zero values.
	 * 
	 * @return the list of non-zero values
	 */
	public double[] values() {
		return data;
	}
	
	/**
	 * Computes the {@code n}-norm of the vector.
	 * 
	 * @param n	 the index of the norm
	 * @return the {@code n}-norm of the vector
	 */
	public double norm(int n) {
		assert(n > 0);
		
		double result = 0;
		for (int i : ids) {
			result += Math.abs(Math.pow(data[i], n));
		}
		result = Math.pow(result, 1.0 / (double) n);
		return result;
	}
	
	@Override
	public SparseVector clone() {
		return new SparseVector(this);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		
		if (!(obj instanceof SparseVector)) {
			return false;
		}
		
		SparseVector that = (SparseVector) obj;
		
		if (size != that.size) {
			return false;
		}
		
		boolean equal = true;
		for (int i : ids) {
			equal &= Misc.isEqual(data[i], that.get(i));
		}
		
		return equal;
	}
	
	@Override
	public int hashCode() {
		return 7 * size + 17 * ids.hashCode() + 31 * data.hashCode();
	}
	
	/**
	 * Converts the sparse vector to an array of doubles.
	 * 
	 * @return the corresponding array of doubles
	 */
	public double[] toArray() {
		double[] result = new double[size];
		for (int i = 0; i < count; ++i) {
			result[ids[i]] = data[i];
		}
		return result;
	}
	
	/**
	 * Converts the sparse vector to a dense representation.
	 * 
	 * @return the dense representation of the vector
	 */
	public DenseVector toDense() {
		DenseVector result = new DenseVector(size);
		for (int i : ids) {
			result.set(i, data[i]);
		}
		return result;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for (int i = 0; i < size; i++) {
			int j = Arrays.binarySearch(ids, i);
			sb.append(j >= 0 ? data[j] : 0.0);
			if (i < size - 1) {
				sb.append(", ");
			}
		}
		sb.append("]");

		return sb.toString();
	}
}
