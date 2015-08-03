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
package aljebra.utils;

/**
 * Class which contains some utility functions, such as functions for sorting,
 * calculating max value in an array and so on.
 * 
 * @author Mirko Polato
 *
 */
public class Misc {

	// Hides the default constructor.
	private Misc() {}

	/**
	 * Returns the max value in the given array of integers.
	 * 
	 * @param array	 the array
	 * @return the maximum value inside the array
	 */
	public static int max(int[] array) {
		assert(array.length > 0);

		int m = array[0];
		for (int i = 1; i < array.length; ++i) {
			if (array[i] > m) {
				m = array[i];
			}
		}
		return m;
	}

	/**
	 * Sorts (using Quicksort) the {@code inputArray} and consequently reorders the {@code compInt} and {@code compDbl}
	 * arrays (if not null).
	 * 
	 * @param start	 	the first index of the portion to sort
	 * @param end 	 	the last index of the portion to sort
	 * @param inputArr	the array to sort
	 * @param compInt	the companion array of integers
	 * @param compDbl	the companion array of doubles
	 */
	public static void sort(int start, int end, int[] inputArr, int[] compInt, double[] compDbl) {
		new QuickSort().sort(start, end, inputArr, compInt, compDbl);
	}

	/**
	 * Sorts (using Quicksort) the {@code inputArray} and consequently reorders the {@code compInt} and {@code compDbl}
	 * arrays (if not null).
	 * 
	 * @param inputArr	the array to sort
	 * @param compInt	the companion array of integers
	 * @param compDbl	the companion array of doubles
	 */
	public static void sort(int[] inputArr, int[] compInt, double[] compDbl) {
		sort(0, inputArr.length - 1, inputArr, compInt, compDbl);
	}

	/**
	 * Returns true whether {@code d1} is in the neighborhood, of radius
	 * 1e-6 of {@code d2}. 
	 * 
	 * @param d1	the first double
	 * @param d2	the second double
	 * @return true whether {@code d1} is in the neighborhood of {@code d2} 
	 */
	public static boolean isEqual(double d1, double d2) {
		return Math.abs(d1 - d2) < 1e-6;
	}

	/**
	 * Calculates the hypotenuse of a right-angle triangle 
	 * from (0, 0) to ({@code x}, {@code y}). 
	 * 
	 * @param x	 the x-axis coordinate 
	 * @param y	 the y-axis coordinate
	 * @return the hypotenuse magnitude.
	 */
	public static double hypot(double x, double y) {
		double r;
		if (Math.abs(x) > Math.abs(y)) {
			r = y / x;
			r = Math.abs(x) * Math.sqrt(1 + r * r);
		} else if (!isEqual(y, 0.0)) {
			r = x / y;
			r = Math.abs(y) * Math.sqrt(1 + r * r);
		} else {
			r = 0.0;
		}
		return r;
	}

	/**
	 * QuickSort support class.
	 */
	private static class QuickSort {

		private int[] array;
		private int[] companionInt;
		private double[] companionDbl;

		public void sort(int start, int end, int[] inputArr, int[] compInt, double[] compDbl) {
			assert(start >= 0 && start < inputArr.length && end < inputArr.length); 

			if (inputArr == null || inputArr.length == 0) {
				return;
			}

			this.array = inputArr;
			this.companionInt = compInt;
			this.companionDbl = compDbl;

			quickSort(start, end);
		}

		private void quickSort(int lowerIndex, int higherIndex) {

			int i = lowerIndex;
			int j = higherIndex;
			int pivot = array[lowerIndex + (higherIndex - lowerIndex) / 2];

			while (i <= j) {
				while (array[i] < pivot) ++i;
				while (array[j] > pivot) --j;

				if (i <= j) {
					exchangeNumbers(i, j);
					++i;
					--j;
				}
			}

			if (lowerIndex < j) {
				quickSort(lowerIndex, j);
			}

			if (i < higherIndex) {
				quickSort(i, higherIndex);
			}
		}

		private void exchangeNumbers(int i, int j) {
			int temp = array[i];
			array[i] = array[j];
			array[j] = temp;

			if (companionInt != null) {
				temp = companionInt[i];
				companionInt[i] = companionInt[j];
				companionInt[j] = temp;
			}

			if (companionDbl != null) {
				double tempd = companionDbl[i];
				companionDbl[i] = companionDbl[j];
				companionDbl[j] = tempd;
			}
		}
	}
}
