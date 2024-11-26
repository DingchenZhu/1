#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include "Circuit.h"

using namespace std;
using namespace Eigen;

bool solve(
    MatrixXd& A,                  // 系统雅可比矩阵
    VectorXd& b,                  // 系统右侧向量
    VectorXd& x,                  // 解向量
    const vector<MOSAdjoint>& mosAdjointList, // MOS伴随模型数据
    int maxIterations = 100,      // 最大迭代次数
    double tolerance = 1e-6       // 收敛阈值
) {
    int nodeNum = A.rows();
    x = VectorXd::Zero(nodeNum);  // 初始解（零向量）

    for (int iter = 0; iter < maxIterations; ++iter) {
        // 更新雅可比矩阵 A 和右侧向量 b
        MatrixXd J = A; // 初始化雅可比矩阵为导纳矩阵
        VectorXd F = b; // 初始化非线性方程值为 b

        for (const auto& mos : mosAdjointList) {
            int pNumD = stoi(mos.drain) + 1;  // 漏极节点编号
            int pNumS = stoi(mos.source) + 1; // 源极节点编号
            int pNumG = stoi(mos.gate) + 1;   // 栅极节点编号

            // 计算非线性贡献
            double vgs = x[pNumG] - x[pNumS];
            double vds = x[pNumD] - x[pNumS];

            // 计算 gm 和 gds
            double gm = mos.gmk;  // 跨导
            double gds = mos.gdsk;  // 漏源导纳
            double ids = gm * vgs + gds * vds + mos.Ik; // 漏源电流

            // 更新雅可比矩阵
            J(pNumD, pNumD) += gds;
            J(pNumD, pNumS) -= gds + gm;
            J(pNumD, pNumG) += gm;

            J(pNumS, pNumD) -= gds;
            J(pNumS, pNumS) += gds + gm;
            J(pNumS, pNumG) -= gm;

            // 更新右侧向量
            F[pNumD] -= ids;  // 漏极电流流出
            F[pNumS] += ids;  // 源极电流流入
        }

        // 线性求解：LU 分解
        VectorXd deltaX = J.lu().solve(-F);

        // 更新解
        x += deltaX;

        // 检查收敛条件
        if (deltaX.norm() < tolerance) {
            cout << "Converged in " << iter + 1 << " iterations." << endl;
            return true;
        }
    }

    cout << "Failed to converge within " << maxIterations << " iterations." << endl;
    return false;
}

