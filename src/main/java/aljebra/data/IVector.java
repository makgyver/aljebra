package aljebra.data;

import java.io.Serializable;

/**
 * Interface for vectors.
 * 
 * @author Mirko Polato
 *
 */
public interface IVector extends Cloneable, Serializable {
	
	/**
	 * Returns the value in position {@code index}.
	 * 
	 * @param index	the index of the entry
	 * @return the value in position {@code index}
	 */
	public double get(int index);
	
	/**
	 * Sets the value in position {@code index} with the 
	 * new {@code value}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the new value
	 */
	public void set(int index, double value);
	
	/**
	 * Sets the entire vector to the given {@code value}.
	 * 
	 * @param value	 the value
	 */
	public void setAll(double value);
	
	/**
	 * Returns the size of the vector.
	 * 
	 * @return the size of the vector
	 */
	public int size();
	
	/**
	 * Returns the opposite vector.
	 * 
	 * @return the opposite vector
	 */
	public IVector opposite();
	
	/**
	 * <p>Addition between vectors.</p>
	 * Adds the vector {@code that}.
	 * 
	 * @param that	the vector to add
	 * @return the resulting vector
	 */
	public IVector add(IVector that);
	
	/**
	 * Adds the given {@code value} to the entry in position
	 * {@code index}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the value to add
	 */
	public void add(int index, double value);
	
	/**
	 * Adds the given {@code value} to all the vector's entries.
	 * 
	 * @param value	 the value to add
	 * @return the resulting vector
	 */
	public IVector add(double value);
	
	/**
	 * <p>Subtraction between vectors.</p>
	 * Subtracts the vector {@code that}.
	 * 
	 * @param that	the vector to subtract
	 * @return the resulting vector
	 */
	public IVector sub(IVector that);
	
	/**
	 * Subtracts the given {@code value} to the entry in position
	 * {@code index}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the value to subtract
	 */
	public void sub(int index, double value);
	
	/**
	 * Subtract the given {@code value} to all the vector's entries.
	 * 
	 * @param value	 the value to subtract
	 * @return the resulting vector
	 */
	public IVector sub(double value);
	
	/**
	 * Scales the vector by the given {@code factor}.
	 * 
	 * @param factor	the scale factor
	 * @return the scaled vector
	 */
	public IVector scale(double factor);
	
	/**
	 * Multiplies the given {@code value} to the entry in position
	 * {@code index}.
	 * 
	 * @param index	 the index of the entry
	 * @param value	 the multiplicative factor
	 */
	public void times(int index, double value);
	
	/**
	 * <p>Point-wise multiplication of vectors.</p>
	 * Multiplies point-wise this vector by the vector {@code that}.
	 * 
	 * @param that	the vector to multiply point-wise
	 * @return the resulting vector
	 */
	public IVector times(IVector that);
	
	/**
	 * <p>Dot product of vectors.</p>
	 * Applies the dot product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting vector
	 */
	public double dot(IVector that);
	
	/**
	 * <p>Outer product of vectors.</p>
	 * Applies the outer product between this and {@code that} vector.
	 * 
	 * @param that	the vector
	 * @return the resulting matrix
	 */
	public IMatrix outer(IVector that);
	
	/**
	 * <p>Product between vector and matrix.</p>
	 * Applies the product between this vector and the given {@code matrix}. </br>
	 * The product between <tt>1 x n</tt> matrix (i.e., n-dimensional vector), 
	 * and an <tt>n x m</tt> matrix is a <tt>1 x m</tt> matrix (i.e., row vector).
	 * 
	 * <p>See also {@link IMatrix#times(IVector)}.</p>
	 * 
	 * @param matrix	the matrix
	 * @return the resulting vector
	 */
	public IVector times(IMatrix matrix);
	
	/**
	 * Computes the sum of the vector's entries.
	 * 
	 * @return the sum of the vector's entries
	 */
	public double sum();
	
	/**
	 * Computes the mean of the vector's entries.
	 * 
	 * @return the mean of the vector's entries
	 */
	public double mean();
	
	/**
	 * Computes the {@code n}-norm of the vector.
	 * 
	 * @param n	 the index of the norm
	 * @return the {@code n}-norm of the vector
	 */
	public double norm(int n);
	
	/**
	 * Converts to array of doubles.
	 * 
	 * @return the corresponding array
	 */
	public double[] toArray();
}
