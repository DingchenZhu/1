#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "Circuit.h"

using namespace std;



void NetListParser(const string& filename, CircuitData& circuitData)
{
    ifstream infile(filename);
    if (!infile.is_open())
    {
        cerr << "Cannot open the file: " << filename << endl;
        return;
    }
    string line;
    while (getline(infile, line))
    {
        if (line.empty() || line[0] == '*')  continue;
        istringstream iss(line);
        string token;
        iss >> token;

        if (token[0] == 'R') { // 电阻
            RINFO r;
            r.Name = token;
            iss >> r.N1 >> r.N2 >> r.Value;
            circuitData.rinfo.push_back(r);
        }
        else if (token[0] == 'C') { // 电容
            CINFO c;
            c.Name = token;
            iss >> c.N1 >> c.N2 >> c.Value;
            circuitData.cinfo.push_back(c);
        }
        else if (token[0] == 'L') { // 电感
            LINFO l;
            l.Name = token;
            iss >> l.N1 >> l.N2 >> l.Value;
            circuitData.linfo.push_back(l);
        }
        else if (token[0] == 'V' || token[0] == 'I') { // 源
            SourceINFO s;
            s.SourceName = token;
            iss >> s.SourceN1 >> s.SourceN2 >> s.Sourcetype;
            if (s.Sourcetype == "DC") {
                iss >> s.SourceDcValue;
            }
            else if (s.Sourcetype == "SIN") {
                iss >> s.SourceDcValue >> s.SourceAcValue >> s.SourceFreq >> s.SourcePhase;
            }
            circuitData.sourceinfo.push_back(s);
        }
        else if (token[0] == 'M') { // MOS管
            MOSINFO m;
            m.MOSName = token;
            iss >> m.MOSN1 >> m.MOSN2 >> m.MOSN3 >> m.MOStype >> m.MOSW >> m.MOSL >> m.MOSID;
            circuitData.mosinfo.push_back(m);
        }
        else if (token == ".MODEL") { // 模型
            ModelINFO model;
            string modelID, keyValue;
            iss >> modelID;
            model.ID = stoi(modelID);
            string str1, str2, str3, str4, str5;
            iss >> str1 >> model.VT >> str2 >> model.MU >> str3 >> model.COX >> str4 >> model.LAMBDA >> str5 >> model.CJ0;
            circuitData.modelinfo.push_back(model);
        }
    }
}

//int main()
//{
//    string filename = "D:\\a_大三上\\模拟EDA\\buffer.sp";
//    CircuitData c1;
//    NetListParser(filename, c1);
//    cout << c1.sourceinfo[0].SourceDcValue << endl;
//}