Aljebra
=======

**Aljebra** is a Java library which contains high performance data structures for linear algebra.


### Code Snippet

This is a simple code example of how to create and use dense and sparse structures.

```java
public void main(String[] args) {

	// Initialize a dense vector and a dense matrix
	DenseVector v = new DenseVector(new double[] {1, 2, 3});
	DenseMatrix m = new DenseMatrix(new double[][] {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});

	// Vector-Matrix multiplication
	DenseVector v_m = v.times(m);

	// Initialize a sparse vector and a sparse matrix
	SparseVector sv = new SparseVector(new int[] {0, 1, 2}, new double[] {1, 2, 3});
	SparseMatrix sm = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2}, 
				new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9});

	// Opposite vector
	SparseVector _sv = sv.opposite();

	// Matrix transposition
	SparseMatrix sm_t = sm.transpose();

}
```


### GPL License

Aljebra is [free software](http://www.gnu.org/philosophy/free-sw.html): you can redistribute it and/or modify it under the terms of the [GNU General Public License (GPL)](http://www.gnu.org/licenses/gpl.html) as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Aljebra is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License along with Aljebra. If not, see http://www.gnu.org/licenses/.
