#include <iostream>
#include "Circuit.h"

int main()
{
	CircuitData circuit;
	string filename = "D:\\a_大三上\\模拟EDA\\buffer.sp";
	NetListParser(filename, circuit);

	circuit = processCircuit(circuit);
	// Test
	cout<<circuit.cinfo[0].N1<<circuit.cinfo[0].N2<<endl;

	// vector<MOSAdjoint> mos_adjoints;
	// mos_adjoints.push_back(ComputeAdjointModel(circuit.mosinfo[0], circuit.modelinfo[0], 2*circuit.sourceinfo[0].SourceDcValue/3, circuit.sourceinfo[0].SourceDcValue/2));
	// mos_adjoints.push_back(ComputeAdjointModel(circuit.mosinfo[1], circuit.modelinfo[1],  2 * circuit.sourceinfo[0].SourceDcValue / 3, circuit.sourceinfo[0].SourceDcValue / 2));
	// mos_adjoints.push_back(ComputeAdjointModel(circuit.mosinfo[2], circuit.modelinfo[0], 2 * circuit.sourceinfo[0].SourceDcValue / 3, circuit.sourceinfo[0].SourceDcValue / 2));
	// mos_adjoints.push_back(ComputeAdjointModel(circuit.mosinfo[3], circuit.modelinfo[1],  2 * circuit.sourceinfo[0].SourceDcValue / 3, circuit.sourceinfo[0].SourceDcValue / 2));

	// MatrixXd A;
	// VectorXd b;
	// vector<string> x;
	// GenMatrix(circuit, 7, A, b, x);

	// VectorXd x_x;
	// bool flag = solve(A, b, x_x, mos_adjoints,100, 1e-6);

	// for(int  i = 0; i < 7; i ++)
	// {
	// 	cout << "Node " << i << " : " << x_x[i] << " V " << endl;
	// }

	// return 0;


}