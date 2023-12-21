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

int main(){
    int N = 100;
    double T = 0;
    char wave = 'd';
    std::string path = "ResultData/1st.json";
    double atol = 1e-6;
    double rtol = 1e-3;
    int maxstep = 1000;
    tJSBMF MFSolver(0, 0, T, N, wave, maxstep, atol, rtol);
    std::vector<double> x_vec = vector_range(0.01, 0.5, 0.01);
    std::vector<double> J_vec = vector_range(0.01, 0.5, 0.01);
    std::vector<std::vector<double> > Delta_vec;
    for(int i=0; i<x_vec.size(); i++){
        std::vector<double> Delta_row;
        for(int j=0; j<J_vec.size(); j++){
            MFSolver.x = x_vec.at(i);
            MFSolver.J = J_vec.at(j);
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
    jsondata["T"] = T;
    jsondata["wave"] = std::string(1, wave);
    jsondata["atol"] = atol;
    jsondata["rtol"] = rtol;
    Json::Value xJson, JJson, DeltaJson, DeltarowJson;
    for(int i=0; i<x_vec.size(); i++){
        xJson.append(x_vec.at(i));
        for(int j=0; j<J_vec.size(); j++){
            JJson.append(J_vec.at(j));
            DeltarowJson.append(Delta_vec.at(i).at(j));
        }
        DeltaJson.append(DeltarowJson);
    }
    jsondata["xList"] = xJson;
    jsondata["JList"] = JJson;
    jsondata["DeltaList"] = DeltaJson;

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