#ifndef _CIRCUIT_H_
#define _CIRCUIT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


struct RINFO {
    string Name;
    string N1;
    string N2;
    double Value;
};

struct CINFO {
    string Name;
    string N1;
    string N2;
    double Value;
};

struct LINFO {
    string Name;
    string N1;
    string N2;
    double Value;
};

struct SourceINFO {
    string SourceName;
    string SourceN1;
    string SourceN2;
    string Sourcetype;
    double SourceDcValue;
    double SourceAcValue;
    double SourceFreq;
    double SourcePhase;
};

struct ControlledSource
{
    string ControlledName;
    string N1;
    string N2;
    string ctrlN1;
    string ctrlN2;
    double value;
    char type;      // E:VCVS,G:VCCS,F:CCCS,H:CCVS
};

struct MOSINFO {
    string MOSName;
    string MOSN1;
    string MOSN2;
    string MOSN3;
    string MOStype;
    double MOSW;
    double MOSL;
    int MOSID;
    string MOSMODEL;
};

struct ModelINFO {
    int ID;
    double VT;
    double MU;
    double COX;
    double LAMBDA;
    double CJ0;
};

// 数据存储结构
struct CircuitData {
    vector<RINFO> rinfo;
    vector<CINFO> cinfo;
    vector<LINFO> linfo;
    vector<SourceINFO> sourceinfo;
    vector<MOSINFO> mosinfo;
    vector<ModelINFO> modelinfo;
};

// 用于伴随模型的结构体
struct MOSAdjoint
{
    string name;
    string drain;
    string source;
    string gate;
    double vds;   // 当前漏-源电压
    double vgs;   // 当前栅-源电压
    double gmk;   // 栅-源跨导
    double gdsk;  // 漏-源导纳
    double Ik;    // 等效电流
};

void calculation_MOS(const MOSINFO& mos, const ModelINFO& model, double Vgs, double Vds, double& gm, double& gds);
void GenMatrix(
    const CircuitData& circuitdata,
    int nodeNum,
    MatrixXd& A,
    VectorXd& b,
    vector<string>& x
);
vector<MOSAdjoint> ProcessMOSFETs(
    const vector<MOSINFO>& mosList,
    const vector<ModelINFO>& modelList,
    const unordered_map<string, pair<double, double>>& voltageMap);
void NetListParser(const string& filename, CircuitData& circuitData);
MOSAdjoint ComputeAdjointModel(const MOSINFO& mos, const ModelINFO& model, double vgs, double vds);
bool solve(
    MatrixXd& A,                  // 系统雅可比矩阵
    VectorXd& b,                  // 系统右侧向量
    VectorXd& x,                  // 解向量
    const vector<MOSAdjoint>& mosAdjointList, // MOS伴随模型数据
    int maxIterations,     // 最大迭代次数
    double tolerance     // 收敛阈值
);
CircuitData processCircuit(CircuitData& circuit) ；
#endif
