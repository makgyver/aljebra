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
package aljebra.data.sparse;

import java.util.ArrayList;
import java.util.Arrays;

import aljebra.data.IMatrix;
import aljebra.data.IVector;
import aljebra.data.dense.DenseVector;
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
public class SparseVector implements IVector {

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
		clean(ids, vals);
		return new SparseVector(ids, vals);
	}
	
	// Clean up the vectors: remove the entry corresponding to zero values.
	private static void clean(int[] ids, double[] v) {
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

		System.arraycopy(vec.data, 0, data, 0, size);
		System.arraycopy(vec.ids, 0, ids, 0, size);
		
		this.count = vec.count;
	}
	
	@Override
	public double get(int index) {
		assert(index >= 0 && index < size);
		
		int i = Arrays.binarySearch(ids, 0, count, index);
		return i >= 0 ? data[i] : 0;
	}
	
	@Override
	public void set(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = getIndex(index);
			data[i] = value;
		}
	}
	
	@Override
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
	
	@Override
	public int size() {
		return size;
	}
	
	@Override
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
	
	@Override
	public SparseVector opposite() {
		double[] newData = new double[count];
		for (int i = 0; i < count; ++i) {
			newData[i] = -data[i];
		}
		return new SparseVector(size, ids, newData);
	}
	
	@Override
	public SparseVector add(IVector that) {
		assert(size == that.size());
		
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) + that.get(i));
		}
		return result;
	}
	
	
	@Override
	public SparseVector add(double value) {
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) + value);
		}
		return result;
	}
	
	@Override
	public void add(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = getIndex(index);
			data[i] += value;
		}
	}
	
	@Override
	public SparseVector sub(IVector that) {
		assert(size == that.size());
		
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			result.set(i, get(i) - that.get(i));
		}
		return result;
	}
	
	@Override
	public SparseVector sub(double value) {
		SparseVector result = new SparseVector(size);
		for (int i : ids) {
			result.set(i, data[i] - value);
		}
		return result;
	}
	
	@Override
	public void sub(int index, double value) {
		assert(index >= 0 && index < size);
		
		if (!Misc.isEqual(value, 0.0)) {
			int i = getIndex(index);
			data[i] -= value;
		}
	}
	
	@Override
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
	
	@Override
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
	
	@Override
	public SparseVector times(IVector that) {
		assert(size == that.size());
		
		SparseVector result = new SparseVector(size);
		for (int i : ids) {
			result.set(i, data[i] * that.get(i));
		}
		return result;
	}
	
	@Override
	public SparseVector div(IVector that) {
		assert(size == that.size());
		
		SparseVector result = new SparseVector(size);
		for (int i : ids) {
			result.set(i, data[i] / that.get(i));
		}
		return result;
	}
	
	@Override
	public SparseVector times(IMatrix matrix) {
		assert(size == matrix.rows());
		
		SparseVector result = new SparseVector(matrix.cols());
		for (int i = 0; i < matrix.cols(); ++i) {
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
	
	@Override
	public double dot(IVector that) {
		assert(size == that.size());
		
		double result = 0.0;
		for (int i : ids) {
			result += data[i] * that.get(i);
		}
		return result;
	}
	
	@Override
	public SparseMatrix outer(IVector that) {
		int[] rows = new int[size * that.size()];
		int[] cols = new int[size * that.size()];
		double[] vals = new double[size * that.size()];
		
		int k = 0;
		for (int i : ids) {
			for (int j = 0; j < that.size(); ++j) {
				if (!Misc.isEqual(that.get(j), 0)) {
					rows[k] = i;
					cols[k] = j;
					vals[k] = data[i] * that.get(j);
					++k;
				}
			}
		}
		return new SparseMatrix(size, that.size(), rows, cols, vals);
	}
	
	@Override
	public double sum() {
		double result = 0;
		for (double v : data) {
			result += v;
		}
		return result;
	}
	
	@Override
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
	
	@Override
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
	
	@Override
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
