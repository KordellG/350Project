#pragma once

#include <iostream>
#include <fstream>

#include "Plotter.h"
#include "Matrix.h"
#include "Vector.h"

/*

Template circuit (showing node indices):

					  (2)
			(0)       L1   (1)
		.---------.---UUU---.
		|  ^      |         |   +
	I1 ( ) |   R1 <      R2 <  vout
		|         |         |   -
		'---------+---------'
				 -+-
	[ v0  ]       '
x = [ v1  ]
	[ iL1 ]

vout = v1

Final Circuit:
					  (2)
			(0)       L1          (1)        R3   (3)
		.---------.---UUU---.------------.---VVV---.
		|  ^      |         |  +         |         |+
	I1 ( ) |   R1 <      R2 < vout   C1 ===       ( ) V1  (4)
		|         |         |  -         |         |
		'---------+---------'------------'---------'
				 -+-
	[ v0  ]       '
	[ v1  ]
x = [ iL1 ]
	[ v3  ]
	[ iV1 ]

vout = v1

*/

using namespace std;

int main()
{
	const double tmax = 1.0;
	const double h = 0.001;

	// add new parameters here:
	const double Ra = 0.5;
	const double La = 10.0e-6;
	const double Jm = 0.1;
	const double Bm = 0.005;
	const double Ke = 0.1;
	const double Kt = 0.1;
	const double v = 10;

	// change the dimensions of the system here:
	Matrix<double> G(8, 8);
	Vector<double> b(8);
	Vector<double> x(8);



	// output:
	ofstream fout;
	Plotter plotter("Final Project");
	fout.open("outfile.csv");
	fout << "torque, eb" << endl;
	plotter.SetLabels("eb(Angular Frequency)", "ia (A)");

	// stamp G matrix here:

	G.initialize(0);
	int va = 0;
	int v2 = 1;
	int eb = 2;
	int vn = 3;
	int wr = 4;
	int ia = 5;
	int iw = 6;
	int it = 7;


	G(va, va) += 1 / Ra;
	G(va, v2) += -1 / Ra;
	G(va, eb) += 0;
	G(va, vn) += 0;
	G(va, wr) += 0;
	G(va, ia) += 0;
	G(va, iw) += 0;
	G(va, it) += 0;


	G(v2, va) += -1 / Ra;
	G(v2, v2) += 1 / Ra;
	G(v2, eb) += 0;
	G(v2, vn) += 0;
	G(v2, wr) += 0;
	G(v2, ia) += 0;
	G(v2, iw) += 0;
	G(v2, it) += 1;


	G(eb, va) += 0;
	G(eb, v2) += 0;
	G(eb, eb) += 0;
	G(eb, vn) += 0;
	G(eb, wr) += 0;
	G(eb, ia) += 1;
	G(eb, iw) += 1;
	G(eb, it) += -1;


	G(vn, va) += 0;
	G(vn, v2) += 0;
	G(vn, eb) += 0;
	G(vn, vn) += 0;
	G(vn, wr) += 0;
	G(vn, ia) += -1;
	G(vn, iw) += -1;
	G(vn, it) += 0;


	G(wr, va) += 0;
	G(wr, v2) += 0;
	G(wr, eb) += 0;
	G(wr, vn) += 0;
	G(wr, wr) += Jm / h + 1 / Bm;
	G(wr, ia) += 0;
	G(wr, iw) += -Kt;
	G(wr, it) += 0;


	G(ia, va) += 0;
	G(ia, v2) += 0;
	G(ia, eb) += 1;
	G(ia, vn) += -1;
	G(ia, wr) += -Ke;
	G(ia, ia) += 0;
	G(ia, iw) += 0;
	G(ia, it) += 0;


	G(iw, va) += 0;
	G(iw, v2) += 0;
	G(iw, eb) += 1;
	G(iw, vn) += -1;
	G(iw, wr) += 0;
	G(iw, ia) += 0;
	G(iw, iw) += 0;
	G(iw, it) += 0;


	G(it, va) += 0;
	G(it, v2) += 1;
	G(it, eb) += -1;
	G(it, vn) += 0;
	G(it, wr) += 0;
	G(it, ia) += 0;
	G(it, iw) += 0;
	G(it, it) += La / h;


	cout << G;
	for (double t = 0.0; t < tmax; t += h)
	{


		plotter.AddRow(t, 0, x[ia]);
		fout << t << "," << x[wr] << endl;
		cout << t << "," << x[wr] << endl;

		// stamp b vector here:

		b.initialize(0); // clear b vector each time step

		b[va] += 1;
		b[v2] += 0;
		b[eb] += 0;
		b[vn] += 0;
		b[wr] += x[wr] * Jm / h;
		b[ia] += 0;
		b[iw] += 0;
		b[it] += -x[it] * La / h;



		// solve the system:
		x = G.computeInverse() * b;
	}

	plotter.Plot();
	fout.close();

	return 0;
}
