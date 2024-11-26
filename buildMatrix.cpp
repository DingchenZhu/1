#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "Circuit.h"

using namespace std;
using namespace Eigen;

//MOSAdjoint ComputeModel(const MOSINFO& mos, const ModelINFO& model, const SourceINFO&source)
MOSAdjoint ComputeAdjointModel(const MOSINFO& mos, const ModelINFO& model, double vgs, double vds) {
	MOSAdjoint adj;
	adj.name = mos.MOSName;
	adj.vds = vds;
	adj.vgs = vgs;
	adj.drain = mos.MOSN1;
	adj.gate = mos.MOSN2;
	adj.source = mos.MOSN3;

	double VT = model.VT;
	double MU = model.MU;
	double COX = model.COX;
	double LAMBDA = model.LAMBDA;
	double W = mos.MOSW;
	double L = mos.MOSL;
	double Kn = MU * COX * (W / L);  // 导通系数

	if (mos.MOStype == "n")
	{
		// cutting
		if (vgs <= VT)
		{
			adj.gmk = 0.0;
			adj.gdsk = 1e-12;
			adj.Ik = 0.0;
		}
		else if (vds < (vgs - VT))
		{
			// 线性区
			adj.gmk = Kn * vds;
			adj.gdsk = Kn * (vgs - VT - vds);
			adj.Ik = Kn * ((vgs - VT) * vds - 0.5 * vds * vds);
		}
		else {
			// 饱和区
			adj.gmk = Kn * (vgs - VT);
			adj.gdsk = Kn * (vds * LAMBDA);
			adj.Ik = 0.5 * Kn * (vgs - VT) * (vgs - VT) * (1 + LAMBDA * vds);
		}
	}
	else if(mos.MOStype == "p")
	{
		double vsg = -vgs;  // PMOS 使用 \( V_{SG} \)
		double vds_p = -vds;  // PMOS 使用 \( V_{SD} \)

		if (vsg <= -VT) {
			// 截止区
			adj.gmk = 0.0;
			adj.gdsk = 1e-12;
			adj.Ik = 0.0;
		}
		else if (vds_p < (vsg + VT)) {
			// 线性区
			adj.gmk = Kn * vds_p;
			adj.gdsk = Kn * (vsg + VT - vds_p);
			adj.Ik = Kn * ((vsg + VT) * vds_p - 0.5 * vds_p * vds_p);
		}
		else {
			// 饱和区
			adj.gmk = Kn * (vsg + VT);
			adj.gdsk = Kn * (vds_p * LAMBDA);
			adj.Ik = 0.5 * Kn * (vsg + VT) * (vsg + VT) * (1 + LAMBDA * vds_p);
		}
	}
	return adj;
}


vector<MOSAdjoint> ProcessMOSFETs(
	const vector<MOSINFO>& mosList, 
	const vector<ModelINFO>& modelList,
	const unordered_map<string, pair<double, double>>& voltageMap)
{
	vector<MOSAdjoint> adjointModels;

	//对每一个mos器件，使用find方法寻找与它ID相同的model
	for (const auto& mos : mosList) {
		auto it = find_if(modelList.begin(), modelList.end(), [&mos](const ModelINFO& model) {
			return model.ID == mos.MOSID;
			});

		if (it != modelList.end()) {
			// 提取电压信息（假设 voltageMap 提供节点电压）
			double vgs = voltageMap.at(mos.MOSN3).first - voltageMap.at(mos.MOSN2).first;  // Vg - Vs
			double vds = voltageMap.at(mos.MOSN1).first - voltageMap.at(mos.MOSN2).first;  // Vd - Vs

			// 计算伴随模型
			MOSAdjoint adj = ComputeAdjointModel(mos, *it, vgs, vds);
			adjointModels.push_back(adj);
		}
	}
	return adjointModels;
}

void GenMatrix(
	const CircuitData& circuitdata,
	const vector<MOSAdjoint>& mosAdjointList,
	int nodeNum,
	MatrixXd& A,
	VectorXd& b,
	vector<string>& x
)
{
	int nodeNums = nodeNum + 1;  // 包含接地节点
	int currentVarIndex = nodeNums;  // 当前额外变量索引
	int INum = 0;

	// 初始化矩阵和向量
	A = MatrixXd::Zero(nodeNums, nodeNums);
	b = VectorXd::Zero(nodeNums);
	x.resize(nodeNums);
	for (int i = 0; i < nodeNums; ++i) {
		x[i] = "v_" + to_string(i);
	}
	//对于电阻
	for (const auto& r : circuitdata.rinfo)
	{
		int pNum1 = stoi(r.N1) + 1;
		int pNum2 = stoi(r.N2) + 1;
		const string& name = r.Name;
		double value = r.Value;
		A(pNum1, pNum1) += 1.0 / value;
		A(pNum1, pNum2) -= 1.0 / value;
		A(pNum2, pNum1) -= 1.0 / value;
		A(pNum2, pNum2) += 1.0 / value;
	}
	for (const auto& source : circuitdata.sourceinfo)
	{
		double svalue = source.SourceDcValue;
		x.push_back("I_" + source.SourceName);
		INum++;
		A.conservativeResize(nodeNums + INum, nodeNums + INum);
		b.conservativeResize(nodeNums + INum);
		A(nodeNums + INum - 1, nodeNums + INum - 1) = 0;
		b(nodeNums + INum - 1) = 0;

		int pNum1 = stoi(source.SourceN1) + 1;
		int pNum2 = stoi(source.SourceN2) + 1;

		A(pNum1, nodeNums + INum - 1) += 1.0;
		A(pNum2, nodeNums + INum - 1) -= 1.0;
		A(nodeNums + INum - 1, pNum1) += 1.0;
		A(nodeNums + INum - 1, pNum2) -= 1.0;
		b(nodeNums + INum - 1) = source.SourceDcValue;
	}
	for (const auto& l : circuitdata.linfo)
	{
		x.push_back("I_" + l.Name);
		INum++;
		A.conservativeResize(nodeNums + INum, nodeNums + INum);
		b.conservativeResize(nodeNums + INum);
		A(nodeNums + INum - 1, nodeNums + INum - 1) = 0;
		b(nodeNums + INum - 1) = 0;

		int pNum1 = stoi(l.N1) + 1;
		int pNum2 = stoi(l.N2) + 1;
		A(pNum1, nodeNums + INum - 1) += 1.0;
		A(pNum2, nodeNums + INum - 1) -= 1.0;
		A(nodeNums + INum - 1, pNum1) += 1.0;
		A(nodeNums + INum - 1, pNum2) -= 1.0;
		b(nodeNums + INum - 1) = 0.0;  // 电压源为 0
	}
	for (const auto& mos : mosAdjointList) {
		int pNumD = stoi(mos.drain) + 1;  // 漏极节点编号
		int pNumS = stoi(mos.source) + 1;  // 源极节点编号
		int pNumG = stoi(mos.gate) + 1;  // 栅极节点编号

		// 独立电流源：直接添加到电流向量中
		b(pNumD) -= mos.Ik;  // 漏极电流流出
		b(pNumS) += mos.Ik;  // 源极电流流入

		// 并联电阻 (Rds)：对应于漏极和源极之间的导纳
		double gds = mos.gdsk;  // 导纳值
		A(pNumD, pNumD) += gds;
		A(pNumS, pNumS) += gds;
		A(pNumD, pNumS) -= gds;
		A(pNumS, pNumD) -= gds;

		// 压控电压源 (gm * Vgs)：通过跨导添加到矩阵中
		double gm = mos.gmk;
		A(pNumD, pNumG) += gm;  // Vgs 影响漏极电流
		A(pNumS, pNumG) -= gm;  // Vgs 影响源极电流
	}
}






//void calculation_MOS(const MOSINFO& mos, const ModelINFO& model, double Vgs, double Vds, double& gm, double& gds)
//{
//	if (mos.MOStype == "n")
//	{
//		if (Vgs <= model.VT)
//		{
//			//截止区
//			gm = 0.0;
//			gds = 0.0;
//		}
//		else if (Vds <= Vgs - model.VT)
//		{
//			//线性区
//			gm = model.MU * model.COX * (mos.MOSW / mos.MOSL) * Vds;
//			gds = model.MU * model.COX * (mos.MOSW / mos.MOSL) * (Vgs - Vds - model.VT);
//		}
//		else
//		{
//			gm = model.MU * model.COX * (mos.MOSW / mos.MOSL) * (Vgs - model.VT) * (1 + model.LAMBDA * Vds);
//			gds = model.LAMBDA * 0.5 * model.MU * model.COX * (mos.MOSW / mos.MOSL) * pow(Vgs - model.VT, 2);
//		}
//	}
//	else if (mos.MOStype == "p")
//	{
//		Vgs = -Vgs;
//		Vds = -Vds;
//		if (Vgs <= -model.VT)
//		{
//			gm = 0.0;
//			gds = 0.0;
//		}
//		else if (Vds <= Vgs + model.VT)
//		{
//			gm = model.MU * model.COX * (mos.MOSW / mos.MOSL) * Vds;
//			gds = model.MU * model.COX * (mos.MOSW / mos.MOSL) * (Vgs - Vds + model.VT);
//		}
//		else
//		{
//			gm = model.MU * model.COX * (mos.MOSW / mos.MOSL) * (Vgs + model.VT) * (1 + model.LAMBDA * Vds);
//			gds = model.LAMBDA * 0.5 * model.MU * model.COX * (mos.MOSW / mos.MOSL) * pow(Vgs + model.VT, 2);
//		}
//	}
//}