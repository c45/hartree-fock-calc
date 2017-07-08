package com;

/**
 * Minimal basis STO-3G calculation on HeH+
 * 
 * Main program which calls HFCALC
 * 
 * A. Szabo and N.S. Ostlund, Modern Quantum Chemistry, 1989.
 */
public class Driver {

	public static void main(String[] args) {
		new Driver();
	}

	public Driver() {
		Common c = new Common();

		// 0 no printing, 1 print converged results, 2 print every iteration
		c.iop = 2;
		// STO-nG calculation (n=1, 2 or 3)
		c.n = 3;
		// bond length (au)
		c.r = 1.4632;
		// slater orbital exponent (function 1)
		c.zeta1 = 2.0925;
		// slater orbital exponent (function 2)
		c.zeta2 = 1.24;
		// atomic number (atom a)
		c.za = 2;
		// atomic number (atom b)
		c.zb = 1;
		
		HartreeFock.calc(c);
	}
}
