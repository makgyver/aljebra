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
package aljebra.unit;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import aljebra.data.sparse.SparseMatrix;
import aljebra.data.sparse.SparseVector;

public class SparseTest {

	protected SparseVector v1;
	protected SparseVector v2;
	
	protected SparseMatrix m1;
	protected SparseMatrix m2;
	protected SparseMatrix m3;
	
	@Before
	public void setUp() throws Exception {
		v1 = new SparseVector(new int[] {0, 1, 2}, new double[] {1, 2, 3});
		v2 = new SparseVector(new int[] {0, 1, 2}, new double[] {4, 5, 6});
		
		m1 = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2}, 
				new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9});
		
		m2 = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {10, 11, 12, 13, 14, 15, 16, 17, 18});
		
		m3 = new SparseMatrix(new int[] {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2},
				new int[] {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
				new double[] {19, 20 , 21, 22, 23, 24, 25, 26, 27, 28, 29, 30});
	}

	@Test
	public void test() {
		SparseVector neg = new SparseVector(new int[] {0, 1, 2}, new double[] {-1, -2, -3});
		SparseVector sum = new SparseVector(new int[] {0, 1, 2}, new double[] {5, 7, 9});
		SparseVector diff = new SparseVector(new int[] {0, 1, 2}, new double[] {-3, -3, -3});
		SparseVector scaled = new SparseVector(new int[] {0, 1, 2}, new double[] {3, 6, 9});
		SparseVector v12 = new SparseVector(new int[] {0, 1, 2}, new double[] {4, 10, 18});
		SparseVector v1m1 = new SparseVector(new int[] {0, 1, 2}, new double[] {30, 36, 42});
		SparseVector m1Diag = new SparseVector(new int[] {0, 1, 2}, new double[] {1, 5, 9});
		SparseVector m1v1 = new SparseVector(new int[] {0, 1, 2}, new double[] {14, 32, 50});
		double vecDot = 32;
		double v1Sum = 6;
		double v1Mean = 2;
		double v1Norm1 = 6;
		double v1Norm2 = Math.sqrt(14);
		
		SparseMatrix mNeg = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2}, 
				new double[] {-1, -2, -3, -4, -5, -6, -7, -8, -9});
		SparseMatrix mSum = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {11, 13, 15, 17, 19, 21, 23, 25, 27});
		SparseMatrix mDiff = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {-9, -9, -9, -9, -9, -9, -9, -9, -9});
		SparseMatrix mDot = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {84, 90, 96, 201, 216, 231, 318, 342, 366});
		SparseMatrix outer = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {4, 5, 6, 8, 10, 12, 12, 15, 18});
		SparseMatrix m1m2 = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {10, 22, 36, 52, 70, 90, 112, 136, 162});
		SparseMatrix m1T = new SparseMatrix(new int[] {0, 0, 0, 1, 1, 1, 2, 2, 2},
				new int[] {0, 1, 2, 0, 1, 2, 0, 1, 2},
				new double[] {1, 4, 7, 2, 5, 8, 3, 6, 9});
		double m1trace = 15;
		double m1norm1 = 45;
		
		assertTrue(neg.equals(v1.opposite()));
		assertTrue(v1.equals(v1.clone()));
		assertTrue(v1.equals(new SparseVector(v1)));
		assertTrue(sum.equals(v1.add(v2)));
		assertTrue(diff.equals(v1.sub(v2)));
		assertEquals(vecDot, v1.dot(v2), 1e-6);
		assertTrue(scaled.equals(v1.scale(3)));
		assertTrue(v12.equals(v1.times(v2)));
		assertTrue(outer.equals(v1.outer(v2)));
		assertTrue(v1m1.equals(v1.times(m1)));
		assertTrue(v1.equals(new SparseVector(v1)));
		assertEquals(v1Sum, v1.sum(), 1e-6);
		assertEquals(v1Mean, v1.mean(), 1e-6);
		assertEquals(v1Norm1, v1.norm(1), 1e-6);
		assertEquals(v1Norm2, v1.norm(2), 1e-6);
		assertEquals(v1.get(0), 1, 1e-6);
		assertEquals(v2.get(2), 6, 1e-6);
		
		assertTrue(mNeg.equals(m1.opposite()));
		assertTrue(m1.equals(m1.clone()));
		assertTrue(m1.equals(new SparseMatrix(m1)));
		assertTrue(m1.isSquare());
		assertTrue(!m3.isSquare());
		assertTrue(m1Diag.equals(m1.diag()));
		assertTrue(mSum.equals(m1.add(m2)));
		assertTrue(mDiff.equals(m1.sub(m2)));
		assertTrue(mDot.equals(m1.dot(m2)));
		assertTrue(m1m2.equals(m1.hadamard(m2)));
		assertTrue(m1v1.equals(m1.times(v1)));
		assertTrue(m1T.equals(m1.transpose()));
		assertEquals(m1.trace(), m1trace, 1e-6);
		assertEquals(m1.norm(1), m1norm1, 1e-6);
		
		v1.add(0, 16);
		assertEquals(v1.get(0), 17, 1e-6);
		
		v1.sub(0, 16);
		assertEquals(v1.get(0), 1, 1e-6);
		
		v1.times(0, 17);
		assertEquals(v1.get(0), 17, 1e-6);
		
		v2.setAll(17);
		assertEquals(v2.get(0), 17, 1e-6);
		assertEquals(v2.get(1), 17, 1e-6);
		assertEquals(v2.get(2), 17, 1e-6);
		
		v2.set(0, 1);
		assertEquals(v2.get(0), 1, 1e-6);
		
		m1.set(0, 0, -2);
		assertEquals(m1.get(0, 0), -2, 1e-6);
	}

}
