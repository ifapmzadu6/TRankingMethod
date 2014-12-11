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
    // RBFの数f
    int rbfCount;
    
    // 入力次元・遅れ時間
    int memoryOfModel;
    
    // RBF(rbfCount, memoryOfModel)
    std::vector<std::vector<double>> rbfs;
    
    // 係数alpha
    std::vector<double> alpha;
    
    
public:
    
    // コンストラクタ
    TermRankingMethod(int rbfCount, int memoryOfModel) : rbfCount(rbfCount), memoryOfModel(memoryOfModel) {
        std::random_device random_device;
        std::mt19937 mt(random_device());
        std::uniform_real_distribution<double> score(0.0, 1.0);
        
        for (int i=0; i<rbfCount; i++) {
            std::vector<double> tmpVector;
            
            for (int j=0; j<memoryOfModel; j++) {
                double random = score(mt);
                tmpVector.push_back(random);
            }
            
            rbfs.push_back(tmpVector);
        }
    }
    
    
    // 実行
    void calculate(std::vector<double> &inputsignal) {
        std::cout << "[TermRankingMethod] Start Calculate" << std::endl;
        
        // 入力の平均
        double averageForInput = averageInVector(inputsignal);
        double denominator = 0.0;
        auto iiter = inputsignal.begin();
        while (iiter != inputsignal.end()) {
            denominator += pow((*iiter - averageForInput), 2);
            ++iiter;
        }
        
        for (int iter=0; iter<rbfCount; iter++) {
            std::cout << "[TermRankingMethod] Calculate K = " << iter << std::endl;
            
            
            // 線形二乗法による最適化
            std::vector<std::vector<double>> A = getPhis(inputsignal);
            std::vector<double> talpha = solveLeastSquaresMethod(A, inputsignal);
            
            
            // エラー計算
            std::vector<double> squaredErrors;
            for (int k=0; k<rbfCount; k++) {
                std::vector<double> predictionOutputs = predictionFunctionOutputs(inputsignal, talpha, k);
                double squaredError = errorOfOneStepPrediction(inputsignal, predictionOutputs, denominator);
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
                        
                        rbfs[k].swap(rbfs[k+1]);
                    }
                }
            }
        }
        
        // 係数alphaの決定
        std::vector<std::vector<double>> A = getPhis(inputsignal);
        alpha = solveLeastSquaresMethod(A, inputsignal);
        
        std::cout << "[TermRankingMethod] End Calculate" << std::endl;
    }
    
    std::vector<double> output(std::vector<double> &inputsignal, int dataCount) {
        std::vector<double> outputs;
        std::vector<double> inputs;
        for (int i=0; i<memoryOfModel; i++) {
            inputs.push_back(inputsignal[i]);
        }
        
        for (int i=0; i<dataCount; i++) {
            double output = 0.0;
            for (int rbfIndex=0; rbfIndex<rbfCount; rbfIndex++) {
                double squeredNorm = 0.0;
                for (int j=0; j<memoryOfModel; j++) {
                    squeredNorm += pow((inputs[j] - rbfs[rbfIndex][j]), 2);
                }
                double spread = 1.0;
                output += alpha[rbfIndex] * exp(- spread * squeredNorm);
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
    std::vector<double> predictionFunctionOutputs(std::vector<double> &inputsignal, std::vector<double> alpha, int &k) {
        std::vector<double> outputs;
        
        for (int i=0; i<inputsignal.size()-memoryOfModel; i++) {
            double output = 0.0;
            for (int rbfIndex=0; rbfIndex<k; rbfIndex++) {
                // ノルムの計算
                double squeredNorm = 0.0;
                for (int j=0; j<memoryOfModel; j++) {
                    squeredNorm += pow((inputsignal[i+j] - rbfs[rbfIndex][j]), 2);
                }
                // RBFの計算
                double spread = 1.0;
                output += alpha[rbfIndex] * exp(- spread * squeredNorm);
            }
            outputs.push_back(output);
        }
        
        return outputs;
    }
    
    std::vector<std::vector<double>> getPhis(std::vector<double> &inputsignal) {
        std::vector<std::vector<double>> phis;
        
        for (int rbfIndex=0; rbfIndex<rbfCount; rbfIndex++) {
            std::vector<double> vector;
            for (int i=0; i<inputsignal.size()-memoryOfModel; i++) {
                // ノルムの計算
                double squeredNorm = 0.0;
                for (int j=0; j<memoryOfModel; j++) {
                    squeredNorm += pow((inputsignal[i+j] - rbfs[rbfIndex][j]), 2);
                }
                // RBFの計算
                double spread = 1.0;
                double output = exp(- spread * squeredNorm);
                
                vector.push_back(output);
            }
            phis.push_back(vector);
        }
        
        return phis;
    }
    
    // ワンステップ誤差を出力
    double errorOfOneStepPrediction(const std::vector<double> &inputsignal, const std::vector<double> &output, double &denominator) {
        double numerator = 0.0;
        auto iiter = inputsignal.begin();
        auto iiter_end = inputsignal.end();
        auto oiter = output.begin();
        while (iiter != iiter_end) {
            double distance = (*oiter) - (*iiter);
            numerator += pow(distance, 2);
            ++iiter; ++oiter;
        }
        
        return numerator / denominator;
    }
    
    // 方程式を解く
    std::vector<double> solveLeastSquaresMethod(std::vector<std::vector<double>> &vA, std::vector<double> &vy) {
        
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
        
        std::vector<double> vo;
        for (int i=0; i<o.size(); i++) {
            vo.push_back(o[i]);
        }
        
        return vo;
    }
    
    // ベクターの平均
    double averageInVector(std::vector<double> vector) {
        double average = 0.0;
        auto iter = vector.begin();
        auto iter_end = vector.end();
        while (iter != iter_end) {
            average += (*iter);
            ++iter;
        }
        return average/vector.size();
    }
    
    
};

#endif /* defined(__TermRankingMethod__TermRankingMethod__) */
