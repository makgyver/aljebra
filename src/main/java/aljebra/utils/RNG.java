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

import java.util.Random;

/**
 * Random Number Generator class.
 * 
 * @author Mirko Polato
 *
 */
public class RNG {
	
	// Hides the default constructor.
	private RNG() {}
	
	/**
	 * Base random number generator object.
	 */
	private static Random rand = new Random(System.currentTimeMillis());

	/**
	 * Sets a new seed for the random object.
	 * 
	 * @param seed	the new random seed
	 */
	public static void seed(long seed) {
		rand = new Random(seed);
	}
	
	/**
	 * Returns a new random integer from a uniform distribution in the range [0, max-1].
	 * 
	 * @param max	the (not included) upper bound of the range 
	 * @return the new random integer number
	 */
	public static int uniformInt(int max) {
		return uniformInt(0, max);
	}

	/**
	 * Returns a new random integer from a uniform distribution in the range [min, max-1].
	 * 
	 * @param min	the (included) lower bound of the range
	 * @param max	the (not included) upper bound of the range 
	 * @return the new random integer number
	 */
	public static int uniformInt(int min, int max) {
		return min + rand.nextInt(max - min);
	}
	
	/**
	 * Returns a new random double from a uniform distribution in the range [0, 1).
	 * 
	 * @return the new random double number
	 */
	public static double uniformDbl() {
		return uniformDbl(0.0, 1.0);
	}

	/**
	 * Returns a new random double from a uniform distribution in the range [min, max-1].
	 * 
	 * @param min	the (included) lower bound of the range
	 * @param max	the (not included) upper bound of the range 
	 * @return the new random double number
	 */
	public static double uniformDbl(double min, double max) {
		return min + (max - min) * rand.nextDouble();
	}
	
	/**
	 * Returns a Bernoulli trial result (true: success, false: failure).
	 * 
	 * @param p	 the probability of success
	 * @return the Bernoulli trial result
	 */
	public static boolean bernoulli(double p) {
		return uniformDbl() < p;
	}
}
