#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include "Circuit.h"
using namespace std;

// 假设已有结构定义


CircuitData processCircuit(CircuitData& circuit) {
    // 1. 提取所有节点名
    set<string> nodes;
    nodes.insert("0"); // 地节点固定为 "0"
    
    for (const auto& r : circuit.rinfo) {
        nodes.insert(r.N1);
        nodes.insert(r.N2);
    }
    for (const auto& c : circuit.cinfo) {
        nodes.insert(c.N1);
        nodes.insert(c.N2);
    }
    for (const auto& l : circuit.linfo) {
        nodes.insert(l.N1);
        nodes.insert(l.N2);
    }
    for (const auto& s : circuit.sourceinfo) {
        nodes.insert(s.SourceN1);
        nodes.insert(s.SourceN2);
    }
    for (const auto& m : circuit.mosinfo) {
        nodes.insert(m.MOSN1);
        nodes.insert(m.MOSN2);
        nodes.insert(m.MOSN3);
    }
    
    // 2. 创建映射表
    map<string, int> nodeMap;
    int nextNodeId = 1;
    for (const auto& node : nodes) {
        if (node == "0") {
            nodeMap[node] = 0; // 地节点固定为 0
        } else {
            nodeMap[node] = nextNodeId++;
        }
    }
    
    for (auto r : circuit.rinfo) {
        r.N1 = to_string(nodeMap[r.N1]);
        r.N2 = to_string(nodeMap[r.N2]);
    }
    for (auto c : circuit.cinfo) {
        c.N1 = to_string(nodeMap[c.N1]);
        c.N2 = to_string(nodeMap[c.N2]);
    }
    for (auto l : circuit.linfo) {
        l.N1 = to_string(nodeMap[l.N1]);
        l.N2 = to_string(nodeMap[l.N2]);
    }
    for (auto s : circuit.sourceinfo) {
        s.SourceN1 = to_string(nodeMap[s.SourceN1]);
        s.SourceN2 = to_string(nodeMap[s.SourceN2]);
    }
    for (auto m : circuit.mosinfo) {
        m.MOSN1 = to_string(nodeMap[m.MOSN1]);
        m.MOSN2 = to_string(nodeMap[m.MOSN2]);
        m.MOSN3 = to_string(nodeMap[m.MOSN3]);
    }
    
    // ControlledSource 更新（如果有）
    for (auto cs : circuit.controlledSource) {
        cs.N1 = to_string(nodeMap[cs.N1]);
        cs.N2 = to_string(nodeMap[cs.N2]);
        cs.ctrlN1 = to_string(nodeMap[cs.ctrlN1]);
        cs.ctrlN2 = to_string(nodeMap[cs.ctrlN2]);
    }
    
    return circuit;
}
