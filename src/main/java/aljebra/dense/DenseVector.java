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

import aljebra.sparse.SparseVector;
import aljebra.utils.Misc;
import aljebra.utils.RNG;

/**
 * <p>Dense implementation of a vector.</p>
 * 
 * <strong>This implementation is a re-elaboration of the one
 * provided by Guibing Guo in his <a href="http://www.librec.net">librec</a>.</strong>
 * 
 * @author Mirko Polato
 *
 */
public class DenseVector implements Cloneable, Serializable {

	private static final long serialVersionUID = -4068089840404586825L;
	
	//
	// STATIC SECTION
	//
	
	/**
	 * Constructs a new dense vector with all entries set to 0.
	 * 
	 * @param size	the size of the vector
	 * @return a new dense null vector
	 */
	public static DenseVector zero(int size) {
		return new DenseVector(size);
	}
	
	/**
	 * Constructs a new dense vector with all entries set to 1.
	 * 
	 * @param size 	the size of the vector
	 * @return a new dense one vector
	 */
	public static DenseVector one(int size) {
		DenseVector one = new DenseVector(size);
		one.setAll(1.0);
		return one;
	}
	
	/**
	 * Constructs a new dense vector initialized with random
	 * double values in the range [0,1].
	 * 
	 * @param size	the size of the vector
	 * @return a new dense vector.
	 */
	public static DenseVector rand(int size) {
		DenseVector vec = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			vec.set(i, RNG.uniformDbl());
		}
		return vec;
	}
	//
	// END STATIC SECTION
	//
	
	
	// The size of the vector
	protected final int size;
	
	// The actual vector's data
	protected double[] data;
	
	/**
	 * Creates a new zero vector.
	 * 
	 * @param size	the size of the vector
	 */
	public DenseVector(int size) {
		this.size = size;
		this.data = new double[size]; 
	}
	
	/**
	 * Creates a new vector with entries taken from the 
	 * given {@code array}.
	 * 
	 * @param array	the data array
	 */
	public DenseVector(double[] array) {
		this.size = array.length;
		this.data = new double[size];
		for (int i = 0; i < size; ++i) {
			this.data[i] = array[i];
		}
	}
	
	/**
	 * Copy constructor.
	 * 
	 * @param vec	the dense vector to copy
	 */
	public DenseVector(DenseVector vec) {
		this(vec.data);
	}
	
	@Override
	public DenseVector clone() {
		return new DenseVector(this);
	}

	/**
	 * Returns the value in position {@code index}.
	 * 
	 * @param index	the index of the entry
	 * @return the value in position {@code index}
	 */
	public double get(int index) {
		assert(index >= 0 && index < size);
		return data[index];
	}
	
	/**
	 * Sets the value in position {@code index} with the 
	 * new {@code value}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the new value
	 */
	public void set(int index, double value) {
		assert(index >= 0 && index < size);
		data[index] = value;
	}

	/**
	 * Sets the entire vector to the given {@code value}.
	 * 
	 * @param value	 the value
	 */
	public void setAll(double value) {
		for (int i = 0; i < size; ++i) {
			data[i] = value;
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
	 * Returns the opposite dense vector.
	 * 
	 * @return the opposite dense vector
	 */
	public DenseVector opposite() {
		DenseVector result = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			result.data[i] = -data[i];
		}
		return result;
	}
	
	/**
	 * <p>Addition between dense vectors.</p>
	 * Adds the vector {@code that}.
	 * 
	 * @param that	the vector to add
	 * @return the resulting vector
	 */
	public DenseVector add(DenseVector that) {
		assert(size == that.size);
		
		DenseVector result = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			result.data[i] = data[i] + that.data[i];
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
		data[index] += value;
	}

	/**
	 * Adds the given {@code value} to all the vector's entries.
	 * 
	 * @param value	 the value to add
	 * @return the resulting vector
	 */
	public DenseVector add(double value) {
		DenseVector result = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			result.data[i] = data[i] + value;
		}
		return result;
	}
	
	/**
	 * <p>Subtraction between dense vectors.</p>
	 * Subtracts the vector {@code that}.
	 * 
	 * @param that	the vector to subtract
	 * @return the resulting vector
	 */
	public DenseVector sub(DenseVector that) {
		assert(size == that.size);
		
		DenseVector result = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			result.data[i] = data[i] - that.data[i];
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
		data[index] -= value;
	}
	
	/**
	 * Subtract the given {@code value} to all the vector's entries.
	 * 
	 * @param value	 the value to subtract
	 * @return the resulting vector
	 */
	public DenseVector sub(double value) {
		DenseVector result = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			result.data[i] = data[i] - value;
		}
		return result;
	}
	
	/**
	 * Scales the vector by the given {@code factor}.
	 * 
	 * @param factor	the scale factor
	 * @return the scaled vector
	 */
	public DenseVector scale(double factor) {
		if (factor == 0.0) {
			return new DenseVector(size);
		}
		
		DenseVector scaled = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			scaled.data[i] = factor * data[i];
		}
		return scaled;
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
		data[index] *= value;
	}
	
	/**
	 * <p>Point-wise multiplication of dense vectors.</p>
	 * Multiplies point-wise this vector by the vector {@code that}.
	 * 
	 * @param that	the vector to multiply point-wise
	 * @return the resulting vector
	 */
	public DenseVector times(DenseVector that) {
		assert(size == that.size);
		
		DenseVector result = new DenseVector(size);
		for (int i = 0; i < size; ++i) {
			result.data[i] = data[i] * that.data[i];
		}
		return result;
	}
	
	/**
	 * <p>Dot product of dense vectors.</p>
	 * Applies the dot product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting vector
	 */
	public double dot(DenseVector that) {
		assert(size == that.size);
		
		double result = 0.0;
		for (int i = 0; i < size; ++i) {
			result += data[i] * that.data[i];
		}
		return result;
	}
	
	/**
	 * <p>Outer product of dense vectors.</p>
	 * Applies the outer product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting dense matrix
	 */
	public DenseMatrix outer(DenseVector that) {
		DenseMatrix result = new DenseMatrix(size, that.size());
		
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < that.size; ++j) {
				result.set(i, j, data[i] * that.data[j]);
			}
		}
		return result;
	}
	
	/**
	 * <p>Product between vector and matrix.</p>
	 * Applies the product between this vector and the given {@code matrix}. </br>
	 * The product between <tt>1 x n</tt> matrix (i.e., n-dimensional vector), 
	 * and an <tt>n x m</tt> matrix is a <tt>1 x m</tt> matrix (i.e., row vector).
	 * 
	 * <p>See also {@link DenseMatrix#times(DenseVector)}.</p>
	 * 
	 * @param matrix	the matrix
	 * @return the resulting vector
	 */
	public DenseVector times(DenseMatrix matrix) {
		assert(size == matrix.rows());
		
		DenseVector result = new DenseVector(matrix.cols());
		for (int i = 0; i < matrix.cols; ++i) {
			double s = 0;
			for (int j = 0; j < size; ++j) {
				s += data[j] * matrix.get(j, i);
			}
			result.data[i] = s;
		}
		return result;
	}
	
	/**
	 * Computes the sum of the vector's entries.
	 * 
	 * @return the sum of the vector's entries
	 */
	public double sum() {
		double result = 0;
		for (int i = 0; i < size; ++i) {
			result += data[i];
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
	 * Converts to array of doubles.
	 * 
	 * @return the corresponding array
	 */
	public double[] toArray() {
		return data.clone();
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
		for (int i = 0; i < size; ++i) {
			result += Math.abs(Math.pow(data[i], n));
		}
		result = Math.pow(result, 1.0 / (double) n);
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		
		if (!(obj instanceof DenseVector)) {
			return false;
		}
		
		DenseVector that = (DenseVector) obj;
		
		if (size != that.size) {
			return false;
		}
		
		boolean equal = true;
		for (int i = 0; i < size && equal; ++i) {
			equal &= Misc.isEqual(data[i], that.data[i]);
		}
		
		return equal;
	}
	
	@Override
	public int hashCode() {
		return 17 * size + data.hashCode();
	}
	
	/**
	 * Converts to sparse vector.
	 * 
	 * @return the corresponding sparse vector
	 */
	public SparseVector toSparse() {
		SparseVector result = new SparseVector(size);
		for (int i = 0; i < size; ++i) {
			if (!Misc.isEqual(data[i], 0.0)) {
				result.set(i, data[i]);
			}
		}
		return result;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for (int i = 0; i < size; ++i) {
			sb.append(data[i]);
			if (i < size - 1)
				sb.append(", ");
		}
		sb.append("]");

		return sb.toString();
	}
}
