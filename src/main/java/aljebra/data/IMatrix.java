package aljebra.data;

import java.io.Serializable;

import aljebra.data.dense.DenseMatrix;
import aljebra.data.dense.DenseVector;

/**
 * Interface for matrices.
 * 
 * @author Mirko Polato
 *
 */
public interface IMatrix extends Cloneable, Serializable {

	/**
	 * Returns the matrix entry in position ({@code row}, {@code col}).
	 * 
	 * @param row	the row of the entry
	 * @param col	the column of the entry
	 * @return	the entry in position {@code row}, {@code col}
	 */
	public double get(int row, int col);
	
	/**
	 * Sets the entry in position ({@code row}, {@code col}) to the given {@code value}.
	 * 
	 * @param row	 the row of the entry
	 * @param col	 the column of the entry
	 * @param value	 the new value for the entry
	 */
	public void set(int row, int col, double value);
	
	//public void setAll(double value);
	
	/**
	 * Returns the number of rows of the matrix.
	 * 
	 * @return the number of rows of the matrix
	 */
	public int rows();
	
	/**
	 * Returns the number of columns of the matrix.
	 * 
	 * @return the number of columns of the matrix
	 */
	public int cols();
	
	/**
	 * Returns true if the matrix is square.
	 * 
	 * @return whether the matrix is square or not
	 */
	public boolean isSquare();
	
	/**
	 * Returns the {@code i}-th row as a vector.
	 * 
	 * @param i	 the row's index
	 * @return the {@code i}-th row as a vector
	 */
	public IVector row(int i);
	
	/**
	 * Returns the {@code j}-th column as a vector.
	 * 
	 * @param j  the column's index
	 * @return the {@code j}-th column as a vector
	 */
	public IVector col(int j);
	
	/**
	 * Returns the main diagonal of the matrix as a vector.
	 * 
	 * @return the main diagonal of the matrix as a vector
	 */
	public IVector diag();
	
	//public void setRow(int row, IVector vec);
	
	//public void setCol(int col, IVector vec);
	
	/**
	 * Returns the opposite matrix.
	 * 
	 * @return the opposite matrix
	 */
	public IMatrix opposite();
	
	/**
	 * Apply the given function to all the entries.
	 * 
	 * @param fun	the function
	 */
	public void apply(IFunction fun);
	
	/**
	 * <p>Addition between matrices.</p>
	 * Adds the matrix {@code that}.
	 * 
	 * @param that	the matrix to add
	 * @return the resulting matrix
	 */
	public IMatrix add(IMatrix that);
	
	/**
	 * Adds the given {@code value} to the entry in position
	 * ({@code row}, {@code col}).
	 * 
	 * @param row	 the entry's row
	 * @param col	 the entry's column
	 * @param value	 the value to add
	 */
	public void add(int row, int col, double value);
	
	/**
	 * Adds the given {@code value} to all the entries.
	 * 
	 * @param value	 the value to add
	 * @return the resulting matrix
	 */
	public IMatrix add(double value);
	
	/**
	 * <p>Subtraction between matrices.</p>
	 * Subtracts the matrix {@code that}.
	 * 
	 * @param that	the matrix to subtract
	 * @return the resulting matrix
	 */
	public IMatrix sub(IMatrix that);
	
	/**
	 * Subtracts the given {@code value} to the entry in position
	 * ({@code row}, {@code col}).
	 * 
	 * @param row	 the entry's row
	 * @param col	 the entry's column
	 * @param value	 the value to subtract
	 */
	public void sub(int row, int col, double value);
	
	/**
	 * Subtracts the given {@code value} to all the entries.
	 * 
	 * @param value	 the value to subtract
	 * @return the resulting matrix
	 */
	public IMatrix sub(double value);
	
	/**
	 * Scales each entry by {@code factor}.
	 * 
	 * @param factor	the scale factor
	 * @return the scaled matrix
	 */
	public IMatrix scale(double factor);
	
	/**
	 * <p>Point-wise multiplication of matrices.</p>
	 * Multiplies point-wise this matrix by the matrix {@code that}.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public IMatrix hadamard(IMatrix that);
	
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
	public IVector times(IVector vec);
	
	/**
	 * <p>Dot product of matrices.</p>
	 * Applies the dot product between this and {@code that} matrix.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public IMatrix dot(IMatrix that);
	
	/**
	 * <p>Point-wise divison of matrices.</p>
	 * Divides point-wise this matrix by the matrix {@code that}.
	 * 
	 * @param that	the matrix
	 * @return the resulting matrix
	 */
	public IMatrix div(IMatrix that);
	
	/**
	 * <p>Point-wise power of matrices.</p>
	 * Exponentiates point-wise this matrix by the given {@code power}.
	 * 
	 * @param power	 the exponent
	 * @return the resulting matrix
	 */
	public IMatrix pow(double power);
	
	/**
	 * Computes the transposition of the matrix.
	 * 
	 * @return the transposition of the matrix
	 */
	public IMatrix transpose();
	
	/**
	 * Computes the trace of the matrix: sum of all values
	 * in the main diagonal.
	 * 
	 * @return the trace of the matrix
	 */
	public double trace();
	
	/**
	 * Computes the sum of all entries.
	 * 
	 * @return the sum of all entries
	 */
	public double sum();
	
	/**
	 * Returns the mean value of the matrix.
	 * 
	 * @return the mean value of the matrix
	 */
	public double mean();
	
	//public IMatrix inv();
	
	/**
	 * Computes the {@code n}-norm of the matrix.
	 * 
	 * @param n		the norm's index
	 * @return the {@code n}-norm of the matrix
	 */
	public double norm(int n);
	
	/**
	 * Converts the matrix to a bi-dimensional vector of doubles.
	 * 
	 * @return the bi-dimensional vector representation of the matrix
	 */
	public double[][] toArray();
	
	Object clone() throws CloneNotSupportedException;
	
	public IMatrix copy();
}
