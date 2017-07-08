package com;

/**
 * Common block
 */
public class Common {

	/**
	 * 0 no printing, 1 print converged results, 2 print every iteration
	 */
	public int iop;

	/**
	 * STO-nG calculation (n=1, 2 or 3)
	 */
	public int n;

	/**
	 * bond length (au)
	 */
	public double r;

	/**
	 * slater orbital exponent (function 1)
	 */
	public double zeta1;

	/**
	 * slater orbital exponent (function 2)
	 */
	public double zeta2;

	/**
	 * atomic number (atom a)
	 */
	public double za;

	/**
	 * atomic number (atom b)
	 */
	public double zb;

	public double s12;

	public double t11;

	public double t12;

	public double t22;

	public double v11a;

	public double v12a;

	public double v22a;

	public double v11b;

	public double v12b;

	public double v22b;

	public double v1111;

	public double v2111;

	public double v2121;

	public double v2211;

	public double v2221;

	public double v2222;

	public double[][] s = new double[2][2];

	public double[][] x = new double[2][2];

	public double[][] xt = new double[2][2];

	public double[][] h = new double[2][2];

	public double[][] f = new double[2][2];

	public double[][] g = new double[2][2];

	public double[][] c = new double[2][2];

	public double[][] fprime = new double[2][2];

	public double[][] cprime = new double[2][2];

	public double[][] p = new double[2][2];

	public double[][] oldp = new double[2][2];

	public double[][][][] tt = new double[2][2][2][2];

	public double[][] e = new double[2][2];
}
