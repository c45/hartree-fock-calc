package com;

/**
 * HartreeFock
 */
public class HartreeFock {

	public static double pi = Math.acos(-1);

	/**
	 * calc
	 * 
	 * @param iop    0 no printing, 1 print converged results, 2 print every iteration
	 * @param n      STO-nG calculation (n=1, 2 or 3)
	 * @param r      bond length (au)
	 * @param zeta1  slater orbital exponent (function 1)
	 * @param zeta2  slater orbital exponent (function 2)
	 * @param za     atomic number (atom a)
	 * @param zb     atomic number (atom b)
	 */
	public static void calc(Common c) {
		System.out.println(String.format("STO-%dG for atomic numbers %5.2f and %5.2f", c.n, c.za, c.zb));
		// compute one and two electron integrals
		intgrl(c);
		collect(c);
		scf(c);
	}

	/**
	 * intgrl
	 * 
	 * @param iop
	 * @param n
	 * @param r
	 * @param zeta1
	 * @param zeta2
	 * @param za
	 * @param zb
	 */
	private static void intgrl(Common c) {
		double[][] coeff = new double[3][3];
		double[][] expon = new double[3][3];
		double[] d1 = new double[3];
		double[] a1 = new double[3];
		double[] d2 = new double[3];
		double[] a2 = new double[3];

		coeff[0][0] = 1.0;
		coeff[1][0] = 0.0;
		coeff[2][0] = 0.0;
		coeff[0][1] = 0.678914;
		coeff[1][1] = 0.430129;
		coeff[2][1] = 0.0;
		coeff[0][2] = 0.444635;
		coeff[1][2] = 0.535328;
		coeff[2][2] = 0.154329;

		expon[0][0] = 0.270950;
		expon[1][0] = 0.0;
		expon[2][0] = 0.0;
		expon[0][1] = 0.151623;
		expon[1][1] = 0.851819;
		expon[2][1] = 0.0;
		expon[0][2] = 0.109818;
		expon[1][2] = 0.405771;
		expon[2][2] = 2.22766;

		double r2 = c.r * c.r;

		for (int i = 0; i < c.n; i++) {
			a1[i] = expon[i][c.n - 1] * c.zeta1 * c.zeta1;
			d1[i] = coeff[i][c.n - 1] * Math.pow(2 * a1[i] / pi, 0.75);
			a2[i] = expon[i][c.n - 1] * c.zeta2 * c.zeta2;
			d2[i] = coeff[i][c.n - 1] * Math.pow(2 * a2[i] / pi, 0.75);
		}

		// calculate one-electron integrals
		for (int i = 0; i < c.n; i++) {
			for (int j = 0; j < c.n; j++) {
				double rap = a2[j] * c.r / (a1[i] + a2[j]);
				double rbp = c.r - rap;
				double rap2 = rap * rap;
				double rbp2 = rbp * rbp;

				c.s12 += s(a1[i], a2[j], r2) * d1[i] * d2[j];

				c.t11 += t(a1[i], a1[j], 0) * d1[i] * d1[j];
				c.t12 += t(a1[i], a2[j], r2) * d1[i] * d2[j];
				c.t22 += t(a2[i], a2[j], 0) * d2[i] * d2[j];

				c.v11a += v(a1[i], a1[j], 0, 0, c.za) * d1[i] * d1[j];
				c.v12a += v(a1[i], a2[j], r2, rap2, c.za) * d1[i] * d2[j];
				c.v22a += v(a2[i], a2[j], 0, r2, c.za) * d2[i] * d2[j];

				c.v11b += v(a1[i], a1[j], 0, r2, c.zb) * d1[i] * d1[j];
				c.v12b += v(a1[i], a2[j], r2, rbp2, c.zb) * d1[i] * d2[j];
				c.v22b += v(a2[i], a2[j], 0, 0, c.zb) * d2[i] * d2[j];
			}
		}

		// calculate two-electron integrals
		for (int i = 0; i < c.n; i++) {
			for (int j = 0; j < c.n; j++) {
				for (int k = 0; k < c.n; k++) {
					for (int l = 0; l < c.n; l++) {
						double rap = a2[i] * c.r / (a2[i] + a1[j]);
						//double rbp = c.r - rap;
						double raq = a2[k] * c.r / (a2[k] + a1[l]);
						double rbq = c.r - raq;
						double rpq = rap - raq;

						double rap2 = rap * rap;
						//double rbp2 = rbp * rbp;
						//double raq2 = raq * raq;
						double rbq2 = rbq * rbq;
						double rpq2 = rpq * rpq;

						c.v1111 += twoe(a1[i], a1[j], a1[k], a1[l], 0, 0, 0) * d1[i] * d1[j] * d1[k] * d1[l];
						c.v2111 += twoe(a2[i], a1[j], a1[k], a1[l], r2, 0, rap2) * d2[i] * d1[j] * d1[k] * d1[l];
						c.v2121 += twoe(a2[i], a1[j], a2[k], a1[l], r2, r2, rpq2) * d2[i] * d1[j] * d2[k] * d1[l];
						c.v2211 += twoe(a2[i], a2[j], a1[k], a1[l], 0, 0, r2) * d2[i] * d2[j] * d1[k] * d1[l];
						c.v2221 += twoe(a2[i], a2[j], a2[k], a1[l], 0, r2, rbq2) * d2[i] * d2[j] * d2[k] * d1[l];
						c.v2222 += twoe(a2[i], a2[j], a2[k], a2[l], 0, 0, 0) * d2[i] * d2[j] * d2[k] * d2[l];
					}
				}
			}
		}

		System.out.println(String.format("%11.6f %11.6f %11.6f %11.6f %11.6f", c.r, c.zeta1, c.zeta2, c.s12, c.t11));
		System.out.println(String.format("%11.6f %11.6f %11.6f %11.6f %11.6f", c.t12, c.t22, c.v11a, c.v12a, c.v22a));
		System.out.println(String.format("%11.6f %11.6f %11.6f %11.6f %11.6f", c.v11b, c.v12b, c.v22b, c.v1111, c.v2111));
		System.out.println(String.format("%11.6f %11.6f %11.6f %11.6f", c.v2121, c.v2211, c.v2221, c.v2222));
	}

	/**
	 * Calculates the f function
	 * f0 only (s-type orbitals)
	 * 
	 * @param arg
	 * @return
	 */
	private static double f0(double arg) {
		double f0;
		if (arg < 1.e-6) {
			f0 = 1 - arg / 3;
		} else {
			f0 = Math.sqrt(pi / arg) * erf(Math.sqrt(arg)) / 2;
		}
		return f0;
	}

	private static double erf(double arg) {
		double p = 0.3275911;
		double[] a = { 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429 };
		double t = 1 / (1 + p * arg);
		double tn = t;
		double poly = a[0] * tn;
		for (int i = 1; i < 5; i++) {
			tn = tn * t;
			poly += a[i] * tn;
		}
		double erf = 1 - poly * Math.exp(-arg * arg);
		return erf;
	}

	/**
	 * Calculates overlaps for un-normalized primitives
	 * 
	 * @param a
	 * @param b
	 * @param rab2
	 * @return
	 */
	private static double s(double a, double b, double rab2) {
		double s = Math.pow(pi / (a + b), 1.5) * Math.exp(-a * b * rab2 / (a + b));
		return s;
	}

	/**
	 * Calculates kinetic energy integrals for un-normalized primitives
	 * 
	 * @param a
	 * @param b
	 * @param rab2
	 * @return
	 */
	private static double t(double a, double b, double rab2) {
		double t = a * b / (a + b) * (3 - 2 * a * b * rab2 / (a + b)) * Math.pow(pi / (a + b), 1.5)
				* Math.exp(-a * b * rab2 / (a + b));
		return t;
	}

	/**
	 * Calculates un-normalized nuclear attraction integrals
	 * 
	 * @param a
	 * @param b
	 * @param rab2
	 * @param rcp2
	 * @param zc
	 * @return
	 */
	private static double v(double a, double b, double rab2, double rcp2, double zc) {
		double v = 2 * pi / (a + b) * f0((a + b) * rcp2) * Math.exp(-a * b * rab2 / (a + b));
		return -v * zc;
	}

	/**
	 * Calculates two-electron integrals for un-normalized primitives
	 * 
	 * @param a
	 * @param b
	 * @param rab2
	 * @param rcp2
	 * @param zc
	 * @return
	 */
	private static double twoe(double a, double b, double c, double d, double rab2, double rcd2, double rpq2) {
		double twoe = 2 * Math.pow(pi, 2.5) / ((a + b) * (c + d) * Math.sqrt(a + b + c + d))
				* f0((a + b) * (c + d) * rpq2 / (a + b + c + d))
				* Math.exp(-a * b * rab2 / (a + b) - c * d * rcd2 / (c + d));
		return twoe;
	}

	/**
	 * Assemble relevant matrices from basic integrals
	 * 
	 * @param a
	 * @param b
	 * @param rab2
	 * @param rcp2
	 * @param zc
	 * @return
	 */
	private static void collect(Common c) {
		// form core Hamiltonian
		c.h[0][0] = c.t11 + c.v11a + c.v11b;
		c.h[0][1] = c.t12 + c.v12a + c.v12b;
		c.h[1][0] = c.h[0][1];
		c.h[1][1] = c.t22 + c.v22a + c.v22b;

		// form overlap matrix
		c.s[0][0] = 1;
		c.s[0][1] = c.s12;
		c.s[1][0] = c.s[0][1];
		c.s[1][1] = 1;

		// use canonical orthogonalization
		c.x[0][0] = 1 / Math.sqrt(2 * (1 + c.s12));
		c.x[1][0] = c.x[0][0];
		c.x[0][1] = 1 / Math.sqrt(2 * (1 - c.s12));
		c.x[1][1] = -c.x[0][1];

		// transpose of transformation matrix
		c.xt[0][0] = c.x[0][0];
		c.xt[0][1] = c.x[1][0];
		c.xt[1][0] = c.x[0][1];
		c.xt[1][1] = c.x[1][1];

		// matrix of two-electron integrals
		c.tt[0][0][0][0] = c.v1111;
		c.tt[1][0][0][0] = c.v2111;
		c.tt[0][1][0][0] = c.v2111;
		c.tt[0][0][1][0] = c.v2111;

		c.tt[0][0][0][1] = c.v2111;
		c.tt[1][0][1][0] = c.v2121;
		c.tt[0][1][1][0] = c.v2121;
		c.tt[1][0][0][1] = c.v2121;

		c.tt[0][1][0][1] = c.v2121;
		c.tt[1][1][0][0] = c.v2211;
		c.tt[0][0][1][1] = c.v2211;
		c.tt[1][1][1][0] = c.v2221;

		c.tt[1][1][0][1] = c.v2221;
		c.tt[1][0][1][1] = c.v2221;
		c.tt[0][1][1][1] = c.v2221;
		c.tt[1][1][1][1] = c.v2222;

		matout(c.s, "S");
		matout(c.x, "X");
		matout(c.h, "H");

		System.out.println();
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					for (int l = 0; l < 2; l++) {
						System.out.println(
								String.format("( %d %d %d %d) %8.6f", i + 1, j + 1, k + 1, l + 1, c.tt[i][j][k][l]));
					}
				}
			}
		}
	}

	/**
	 * 
	 * 
	 * @param c
	 */
	private static void scf(Common c) {

		double crit = 1e-4;
		double delta = 1;
		double en = 0;

		int maxit = 25;
		int iter = 0;

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				c.p[i][j] = 0;
			}
		}

		matout(c.p, "P");

		while (++iter <= maxit && delta > crit) {

			System.out.println();
			System.out.println("START OF ITERATION NUMBER = " + iter);

			formg(c);

			matout(c.g, "G");

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					c.f[i][j] = c.h[i][j] + c.g[i][j];
				}
			}

			en = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					en += c.p[i][j] * (c.h[i][j] + c.f[i][j]) / 2;
				}
			}

			matout(c.f, "F");

			System.out.println();
			System.out.println("ELECTRONIC ENERGY = " + en);

			mult(c.f, c.x, c.g, 2);
			mult(c.xt, c.g, c.fprime, 2);
			diag(c.fprime, c.cprime, c.e);
			mult(c.x, c.cprime, c.c, 2);

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					c.oldp[i][j] = c.p[i][j];
					c.p[i][j] = 0;
					for (int k = 0; k < 1; k++) {
						c.p[i][j] += 2 * c.c[i][k] * c.c[j][k];
					}
				}
			}

			matout(c.fprime, "F'");
			matout(c.cprime, "C'");
			matout(c.e, "E");
			matout(c.c, "C");
			matout(c.p, "P");

			delta = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					delta += (c.p[i][j] - c.oldp[i][j]) * (c.p[i][j] - c.oldp[i][j]);
				}
			}
			delta = Math.sqrt(delta / 4);

			System.out.println();
			System.out.println("DELTA(CONVERGENCE OF DENSITY MATRIX) = " + delta);
		}

		if (iter >= maxit) {
			System.out.println("NO CONVERGENCE IN SCF");
			System.exit(0);
		}

		double ent = en + c.za * c.zb / c.r;

		System.out.println();
		System.out.println("CALCULATION CONVERGED");
		System.out.println();
		System.out.println(String.format("ELECTRONIC ENERGY = %19.11f", en));
		System.out.println();
		System.out.println(String.format("TOTAL ENERGY      = %19.11f", ent));

		mult(c.p, c.s, c.oldp, 2);
		matout(c.oldp, "P");
	}

	/**
	 * Calculates the G matrix from the density matrix and two-electron integrals.
	 * 
	 * @param c
	 */
	private static void formg(Common c) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				c.g[i][j] = 0;
				for (int k = 0; k < 2; k++) {
					for (int l = 0; l < 2; l++) {
						c.g[i][j] += c.p[k][l] * (c.tt[i][j][k][l] - c.tt[i][l][k][j] / 2);
					}
				}
			}
		}
	}

	/**
	 * Diagonalizes f to give eigenvectors in c and e.
	 * theta is the angle describing the solution
	 * 
	 * @param f
	 * @param c
	 * @param e
	 */
	private static void diag(double[][] f, double[][] c, double[][] e) {
		double theta;
		double temp;

		if (Math.abs(f[0][0] - f[1][1]) < 1.e-20) {
			theta = pi / 4;
		} else {
			theta = .5 * Math.atan(2 * f[0][1] / (f[0][0] - f[1][1]));
		}
		c[0][0] = Math.cos(theta);
		c[1][0] = Math.sin(theta);
		c[0][1] = Math.sin(theta);
		c[1][1] = -Math.cos(theta);

		e[0][0] = f[0][0] * Math.cos(theta) * Math.cos(theta) + f[1][1] * Math.sin(theta) * Math.sin(theta)
				+ f[0][1] * Math.sin(2 * theta);
		e[1][1] = f[1][1] * Math.cos(theta) * Math.cos(theta) + f[0][0] * Math.sin(theta) * Math.sin(theta)
				- f[0][1] * Math.sin(2 * theta);
		e[1][0] = 0;
		e[0][1] = 0;

		if (e[1][1] <= e[0][0]) {
			temp = e[1][1];
			e[1][1] = e[0][0];
			e[0][0] = temp;
			temp = c[0][1];
			c[0][1] = c[0][0];
			c[0][0] = temp;
			temp = c[1][1];
			c[1][1] = c[1][0];
			c[1][0] = temp;
		}
	}

	/**
	 * multiplies two matrices
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param m
	 * @return
	 */
	private static void mult(double[][] a, double[][] b, double[][] c, int m) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				c[i][j] = 0;
				for (int k = 0; k < m; k++) {
					c[i][j] += a[i][k] * b[k][j];
				}
			}
		}
	}

	private static void matout(double[][] a, String label) {
		System.out.println();
		System.out.println("THE " + label + " ARRAY");
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a.length; j++) {
				System.out.print(String.format("%9.9e ", a[i][j]));
			}
			System.out.println();
		}
	}
}
