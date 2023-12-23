#include "tJSBMF.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <json/json.h>
#include <string>

template <typename datatype>
int vector_range(std::vector<datatype> *vec, datatype start, datatype stop, datatype step);
template <typename datatype>
std::vector<datatype> vector_range(datatype start, datatype stop, datatype step);
int xJ_phasediagram(double T=0);
int xT_phasediagram(double J=0.3);
int JT_phasediagram(double x=0.1);

int main(){
    xJ_phasediagram();
    // tJSBMF example(0.3, 0, 0, 100, 'd');
    // example.reset();
    // example.self_consistent();
    // std::cout<<example.h<<" "<<example.B<<" "<<example.Delta<<std::endl;
    return 0;
}

int JT_phasediagram(double x){
    int N = 100;
    char wave = 'd';
    std::string path = "ResultData/JT_d.json";
    double atol = 1e-6;
    double rtol = 1e-3;
    int maxstep = 1000;
    tJSBMF MFSolver(0, x, 0, N, wave, maxstep, atol, rtol);
    //test
    // MFSolver.x = 0.41;
    // MFSolver.J = 0.46;
    // MFSolver.reset();
    // MFSolver.self_consistent();
    // std::cout<<MFSolver.x<<" "<<MFSolver.J<<" "<<std::fixed<<std::setprecision(10)<<MFSolver.Delta<<std::endl;

    std::vector<double> T_vec = vector_range(0.002, 0.1, 0.002);
    std::vector<double> J_vec = vector_range(0.01, 0.5, 0.01);
    std::vector<std::vector<double> > Delta_vec;
    for(int i=0; i<J_vec.size(); i++){
        std::vector<double> Delta_row;
        for(int j=0; j<T_vec.size(); j++){
            MFSolver.J = J_vec.at(i);
            MFSolver.T = T_vec.at(j);
            MFSolver.reset();
            MFSolver.self_consistent();
            Delta_row.push_back(MFSolver.Delta);
        }
        Delta_vec.push_back(Delta_row);
    }
    //output the datafile
    std::cout<<"Start to write the json file!"<<std::endl;
    Json::Value jsondata;
    jsondata["N"] = N;
    jsondata["x"] = x;
    jsondata["wave"] = std::string(1, wave);
    jsondata["atol"] = atol;
    jsondata["rtol"] = rtol;
    Json::Value TJson, JJson, DeltaJson;
    for(int i=0; i<J_vec.size(); i++){
        Json::Value DeltaJsonRow;
        JJson.append(J_vec.at(i));
        for(int j=0; j<T_vec.size(); j++){
            if(i==0){
                TJson.append(T_vec.at(j));
            }
            DeltaJsonRow.append(Delta_vec.at(i).at(j));
        }
        DeltaJson.append(DeltaJsonRow);
    }
    jsondata["JList"] = JJson;
    jsondata["TList"] = TJson;
    jsondata["DeltaList"] = DeltaJson;

    std::ofstream jsonfile(path);
    jsonfile<<jsondata;
    return 0;
}
int xT_phasediagram(double J){
    int N = 100;
    char wave = 'd';
    std::string path = "ResultData/xT_d.json";
    double atol = 1e-6;
    double rtol = 1e-3;
    int maxstep = 1000;
    tJSBMF MFSolver(J, 0, 0, N, wave, maxstep, atol, rtol);
    //test
    // MFSolver.x = 0.41;
    // MFSolver.J = 0.46;
    // MFSolver.reset();
    // MFSolver.self_consistent();
    // std::cout<<MFSolver.x<<" "<<MFSolver.J<<" "<<std::fixed<<std::setprecision(10)<<MFSolver.Delta<<std::endl;

    std::vector<double> x_vec = vector_range(0.01, 0.5, 0.01);
    std::vector<double> T_vec = vector_range(0.002, 0.1, 0.002);
    std::vector<std::vector<double> > Delta_vec;
    for(int i=0; i<x_vec.size(); i++){
        std::vector<double> Delta_row;
        for(int j=0; j<T_vec.size(); j++){
            MFSolver.x = x_vec.at(i);
            MFSolver.T = T_vec.at(j);
            MFSolver.reset();
            MFSolver.self_consistent();
            Delta_row.push_back(MFSolver.Delta);
        }
        Delta_vec.push_back(Delta_row);
    }
    //output the datafile
    std::cout<<"Start to write the json file!"<<std::endl;
    Json::Value jsondata;
    jsondata["N"] = N;
    jsondata["J"] = J;
    jsondata["wave"] = std::string(1, wave);
    jsondata["atol"] = atol;
    jsondata["rtol"] = rtol;
    Json::Value xJson, TJson, DeltaJson;
    for(int i=0; i<x_vec.size(); i++){
        Json::Value DeltaJsonRow;
        xJson.append(x_vec.at(i));
        for(int j=0; j<T_vec.size(); j++){
            if(i==0){
                TJson.append(T_vec.at(j));
            }
            DeltaJsonRow.append(Delta_vec.at(i).at(j));
        }
        DeltaJson.append(DeltaJsonRow);
    }
    jsondata["xList"] = xJson;
    jsondata["TList"] = TJson;
    jsondata["DeltaList"] = DeltaJson;

    std::ofstream jsonfile(path);
    jsonfile<<jsondata;
    return 0;
}
int xJ_phasediagram(double T){
    int N = 100;
    char wave = 'd';
    std::string path = "ResultData/xJBd_wave_revised.json";
    double atol = 1e-6;
    double rtol = 1e-3;
    int maxstep = 1000;
    tJSBMF MFSolver(0, 0, T, N, wave, maxstep, atol, rtol);
    //test
    // MFSolver.x = 0.41;
    // MFSolver.J = 0.46;
    // MFSolver.reset();
    // MFSolver.self_consistent();
    // std::cout<<MFSolver.x<<" "<<MFSolver.J<<" "<<std::fixed<<std::setprecision(10)<<MFSolver.Delta<<std::endl;

    std::vector<double> x_vec = vector_range(0., 0.5, 0.01);
    std::vector<double> J_vec = vector_range(0., 0.5, 0.01);
    std::vector<std::vector<double> > Delta_vec, B_vec, DeltaSC_vec;
    for(int i=0; i<x_vec.size(); i++){
        std::vector<double> Delta_row, B_row, DeltaSC_row;
        for(int j=0; j<J_vec.size(); j++){
            MFSolver.x = x_vec.at(i);
            MFSolver.J = J_vec.at(j);
            MFSolver.reset();
            MFSolver.self_consistent();
            Delta_row.push_back(MFSolver.Delta);
            B_row.push_back(MFSolver.B);
            DeltaSC_row.push_back(MFSolver.DeltaSC());
        }
        Delta_vec.push_back(Delta_row);
        B_vec.push_back(B_row);
        DeltaSC_vec.push_back(DeltaSC_row);
    }
    //output the datafile
    std::cout<<"Start to write the json file!"<<std::endl;
    Json::Value jsondata;
    jsondata["N"] = N;
    jsondata["T"] = T;
    jsondata["wave"] = std::string(1, wave);
    jsondata["atol"] = atol;
    jsondata["rtol"] = rtol;
    Json::Value xJson, JJson, DeltaJson, BJson, DeltaSCJson;
    for(int i=0; i<x_vec.size(); i++){
        Json::Value DeltaJsonRow;
        Json::Value BJsonRow, DeltaSCJsonRow;
        xJson.append(x_vec.at(i));
        for(int j=0; j<J_vec.size(); j++){
            if(i==0){
                JJson.append(J_vec.at(j));
            }
            DeltaJsonRow.append(Delta_vec.at(i).at(j));
            BJsonRow.append(B_vec.at(i).at(j));
            DeltaSCJsonRow.append(DeltaSC_vec.at(i).at(j));
        }
        DeltaJson.append(DeltaJsonRow);
        BJson.append(BJsonRow);
        DeltaSCJson.append(DeltaSCJsonRow);
    }
    jsondata["xList"] = xJson;
    jsondata["JList"] = JJson;
    jsondata["DeltaList"] = DeltaJson;
    jsondata["BList"] = BJson;
    jsondata["DeltaSCList"] = DeltaSCJson;

    std::ofstream jsonfile(path);
    jsonfile<<jsondata;
    return 0;
}

template <typename datatype>
int vector_range(std::vector<datatype> *vec, datatype start, datatype stop, datatype step){
    //vec should be void. or the range will be appended to it
    for(datatype elem = start; elem <= stop; elem += step){
        vec->push_back(elem);
    }
    return 0;
}
template <typename datatype>
std::vector<datatype> vector_range(datatype start, datatype stop, datatype step){
    std::vector<datatype> res;
    vector_range(&res, start, stop, step);
    return res;
}