//
//  TermRankingMethod.h
//  TermRankingMethod
//
//  Created by Keisuke Karijuku on 2014/12/01.
//  Copyright (c) 2014年 Keisuke Karijuku. All rights reserved.
//

#ifndef __TermRankingMethod__TermRankingMethod__
#define __TermRankingMethod__TermRankingMethod__

#include <iostream>
#include <vector>
#include <random>

#include "Eigen/Dense"

class TermRankingMethod {
public:

    
    // RBFの数f
    int rbfCount;
    
    // 入力次元・遅れ時間
    int memoryOfModel;
    
    // RBF(rbfCount, memoryOfModel)
    std::vector<std::vector<double>> rbfs;
    
    // RBFの係数
    double spread;
    
    // TermRanking係数alpha
    std::vector<double> alpha;
    
    // 入力データ数
    double dataCount;
    
    
    // コンストラクタ
    TermRankingMethod(int rbfCount, int memoryOfModel) : rbfCount(rbfCount), memoryOfModel(memoryOfModel) {
        std::random_device random_device;
        std::mt19937 mt(random_device());
        std::uniform_real_distribution<double> score(-1.1, 1.1);
        for (int i=0; i<rbfCount; i++) {
            std::vector<double> tmpVector;
            for (int j=0; j<memoryOfModel; j++) {
                double random = score(mt);
                tmpVector.push_back(random);
            }
            rbfs.push_back(tmpVector);
        }
        spread = 1.0;
    }
    
    
    // 実行
    void calculate(const std::vector<double> &inputsignal) {
        std::cout << "[TermRankingMethod] Start Calculate" << std::endl;
        
        dataCount = inputsignal.size();
        
        // 入力の平均
        double averageForInput = averageInVector(inputsignal);
        // 入力の平均と入力の差分
        double diffInputsignalAndAverage = diffBetweenVectorAndValue(inputsignal, averageForInput);
        
        
        for (int iter=0; iter<rbfCount; iter++) {
            std::cout << "[TermRankingMethod] Calculate K = " << iter << std::endl;
            
            
            // 線形二乗法による最適化
            std::vector<std::vector<double>> A = getPhis(inputsignal);
            std::vector<double> talpha = solveLeastSquaresMethod(A, inputsignal);
            
            
            // エラー計算
            std::vector<double> squaredErrors;
            for (int k=0; k<rbfCount; k++) {
                std::vector<double> predictionOutputs = predictionFunctionOutputs(inputsignal, talpha, k);
                double squaredError = errorOfOneStepPrediction(inputsignal, predictionOutputs, diffInputsignalAndAverage);
                squaredErrors.push_back(squaredError);
            }
            
            
            // 差分計算
            std::vector<double> informationDiscrepancys;
            double firstInformationDiscrepancys = -1.0/inputsignal.size();
            informationDiscrepancys.push_back(firstInformationDiscrepancys);
            for (int k=0; k<rbfCount-1; k++) {
                double informationDiscrepancy = log(squaredErrors[k]/squaredErrors[k+1]);
                informationDiscrepancys.push_back(informationDiscrepancy);
            }
            
            
            // 並び替え＆先頭固定
            for (int tk=iter; tk<rbfCount-iter; tk++) {
                for (int k=tk; k<rbfCount-1-tk; k++) {
                    if (informationDiscrepancys[k] < informationDiscrepancys[k+1]) {
                        double tmpA = informationDiscrepancys[k];
                        informationDiscrepancys[k] = informationDiscrepancys[k+1];
                        informationDiscrepancys[k+1] = tmpA;
                        
                        auto tmp = rbfs[k];
                        rbfs[k] = rbfs[k+1];
                        rbfs[k+1] = tmp;
                    }
                }
            }
        }
        
        
        // 最適化終了後、係数alphaの保存
        std::vector<std::vector<double>> A = getPhis(inputsignal);
        alpha = solveLeastSquaresMethod(A, inputsignal);
        
        
        std::cout << "[TermRankingMethod] End Calculate" << std::endl;
    }
    
    std::vector<double> output(const std::vector<double> &inputsignal, const int &dataCount) const {
        std::vector<double> outputs;
        std::vector<double> inputs;
        for (int i=0; i<memoryOfModel; i++) {
            inputs.push_back(inputsignal[i]);
        }
        
        for (int i=0; i<dataCount; i++) {
            double output = 0.0;
            for (int rbfIndex=0; rbfIndex<rbfCount; rbfIndex++) {
                if (rbfIndex==0) {
                    output += alpha[0];
                }
                else {
                    double squeredNorm = 0.0;
                    for (int j=0; j<memoryOfModel; j++) {
                        squeredNorm += pow((inputs[i+j]-rbfs[rbfIndex][j]), 2);
                    }
                    output += alpha[rbfIndex] * exp(-spread*squeredNorm);
                }
            }
            
            for (int t=memoryOfModel-1; t>0; t--) {
                inputs[t] = inputs[t-1];
            }
            inputs[0] = output;
            
            outputs.push_back(output);
        }
        return outputs;
    }
    
    // k番目の予測関数出力
    std::vector<double> predictionFunctionOutputs(const std::vector<double> &inputsignal, const std::vector<double> &alpha, const int &k) const {
        std::vector<double> outputs;
        for (int i=0; i<dataCount-memoryOfModel; i++) {
            double output = 0.0;
            for (int rbfIndex=0; rbfIndex<k; rbfIndex++) {
                if (rbfIndex==0) {
                    output += alpha[0];
                }
                else {
                    double squeredNorm = 0.0;
                    for (int j=0; j<memoryOfModel; j++) {
                        squeredNorm += pow((inputsignal[i+j]-rbfs[rbfIndex][j]), 2);
                    }
                    output += alpha[rbfIndex] * exp(-spread*squeredNorm);
                }
            }
            outputs.push_back(output);
        }
        return outputs;
    }
    
    std::vector<std::vector<double>> getPhis(const std::vector<double> &inputsignal) const {
        std::vector<std::vector<double>> phis;
        for (int rbfIndex=0; rbfIndex<rbfCount; rbfIndex++) {
            std::vector<double> vector;
            for (int i=0; i<dataCount-memoryOfModel; i++) {
                double output;
                if (rbfIndex==0) {
                    output = 1.0;
                }
                else {
                    double squeredNorm = 0.0;
                    for (int j=0; j<memoryOfModel; j++) {
                        squeredNorm += pow((inputsignal[i+j]-rbfs[rbfIndex][j]), 2);
                    }
                    output = exp(-spread*squeredNorm);
                }
                
                vector.push_back(output);
            }
            phis.push_back(vector);
        }
        return phis;
    }
    
    // ワンステップ誤差を出力
    double errorOfOneStepPrediction(const std::vector<double> &inputsignal, const std::vector<double> &output, double &denominator) const {
        double numerator = 0.0;
        auto iiter = inputsignal.begin();
        auto iiter_end = inputsignal.end();
        auto oiter = output.begin();
        while (iiter != iiter_end) {
            numerator += pow((*oiter-*iiter), 2);
            ++iiter; ++oiter;
        }
        return numerator/denominator;
    }
    
    
    std::vector<double> solveLeastSquaresMethod(const std::vector<std::vector<double>> &vA, const std::vector<double> &vy) const {
        Eigen::MatrixXd A(vA[0].size(), vA.size());
        Eigen::VectorXd y(vy.size()-memoryOfModel);
        
        for (int i=0; i<A.cols(); i++) {
            for (int j=0; j<A.rows(); j++) {
                A(j, i) = vA[i][j];
            }
        }
        for (int i=0; i<vy.size()-memoryOfModel; i++) {
            y(i) = vy[i];
        }
        
        Eigen::MatrixXd gtg = A.transpose()*A;
        Eigen::MatrixXd gty = A.transpose()*y;
        // Gt*G*o = Gt*y
        Eigen::VectorXd o = (gtg).colPivHouseholderQr().solve(gty);
        
        std::vector<double> vo(o.size());
        for (int i=0; i<o.size(); i++) {
            vo[i] = o[i];
        }
        return vo;
    }
    
    double averageInVector(const std::vector<double> &vector) const {
        double average = 0.0;
        auto iter = vector.begin();
        auto iter_end = vector.end();
        while (iter != iter_end) {
            average += (*iter);
            ++iter;
        }
        return average/vector.size();
    }
    
    double diffBetweenVectorAndValue(const std::vector<double> &vector, const double &value) const {
        double diff = 0.0;
        auto iter = vector.begin();
        auto iter_end = vector.end();
        while (iter != iter_end) {
            diff += pow((*iter-value), 2);
            ++iter;
        }
        return diff;
    }
    
    // MARK: File Input
    void readAlpha(std::string fileName) {
        std::vector<double> vector;
        std::ifstream alphaIstream(fileName);
        for (int i=0; i<rbfCount; i++) {
            std::string string;
            std::getline(alphaIstream, string);
            double value = std::stof(string);
            vector.push_back(value);
        }
        alpha = vector;
    }
    
    void readRBFs(std::string fileName) {
        std::vector<std::vector<double>> matrix;
        std::ifstream rbfsIstream(fileName);
        for (int i=0; i<rbfCount; i++) {
            std::vector<double> vector;
            for (int j=0; j<memoryOfModel; j++) {
                std::string string;
                std::getline(rbfsIstream, string, ',');
                double value = std::stof(string);
                vector.push_back(value);
            }
            matrix.push_back(vector);
        }
        rbfs = matrix;
    }
    
    // MARK: File output
    
    void writeAlpha(std::string fileName) const {
        std::ofstream alphaFstream(fileName);
        for (int i=0; i<rbfCount; i++) {
            alphaFstream << alpha[i] << std::endl;
        }
        alphaFstream.close();
    }
    
    void writeRBFs(std::string fileName) const {
        std::ofstream rbfsFstream(fileName);
        for (int i=0; i<rbfCount; i++) {
            for (int j=0; j<memoryOfModel; j++) {
                rbfsFstream << rbfs[i][j] << ",";
            }
        }
        rbfsFstream.close();
    }
};

#endif /* defined(__TermRankingMethod__TermRankingMethod__) */
