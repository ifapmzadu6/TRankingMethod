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
#include "KMeansMethod.h"

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
    TermRankingMethod(int rbfCount, int memoryOfModel, double spread) : rbfCount(rbfCount), memoryOfModel(memoryOfModel), spread(spread) {
    }
    
    
    // 実行
    void calculate(const std::vector<double> &inputsignal) {
        
        KMeansMethod kmeans;
        std::vector<std::vector<double>> cluster = kmeans.calculate(inputsignal, memoryOfModel, rbfCount-1);
        rbfs = cluster;
        
        std::cout << "[TermRankingMethod] Start Calculate" << std::endl;
        
        dataCount = inputsignal.size() - memoryOfModel;
        
        // 入力の平均
        double averageForInput = averageInVector(inputsignal);
        // 入力の平均と入力の差分
        double diffInputsignalAndAverage = diffBetweenVectorAndValue(inputsignal, averageForInput);
        
        for (int iter=0; iter<rbfCount; iter++) {
            std::cout << "[TermRankingMethod] Calculate K = " << iter << std::endl;
            
            // 線形二乗法による最適化
            std::vector<std::vector<double>> alphas;
            for (int k=0; k<rbfCount; k++) {
                std::vector<std::vector<double>> A = getPhis(inputsignal, k, dataCount);
                std::vector<double> talpha = solveLeastSquaresMethod(A, inputsignal);
                alphas.push_back(talpha);
            }
            
            
            // エラー計算
            std::vector<double> squaredErrors;
            for (int k=0; k<rbfCount-1; k++) {
                std::vector<double> predictionOutputs = predictionFunctionOutputs(inputsignal, alphas[k], k, dataCount);
                double error = errorOfOneStepPrediction(inputsignal, predictionOutputs, diffInputsignalAndAverage);
                squaredErrors.push_back(error);
            }
            
            
            // 差分計算
            std::vector<double> informationDiscrepancys;
            double firstInformationDiscrepancys = -1.0/inputsignal.size();
            informationDiscrepancys.push_back(firstInformationDiscrepancys);
            for (int k=0; k<rbfCount-2; k++) {
                double informationDiscrepancy = log(squaredErrors[k]/squaredErrors[k+1]);
                informationDiscrepancys.push_back(informationDiscrepancy);
            }
            
            
            // 並び替え＆先頭固定
            for (int tk=iter; tk<rbfCount-iter-2; tk++) {
                for (int k=iter; k<rbfCount-iter-2; k++) {
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
        
        
        std::cout << "[TermRankingMethod] End Calculate" << std::endl;
    }
    
    std::vector<double> output(const std::vector<double> &inputsignal, const int rbfCount, const int &dataCount) {
        std::vector<double> outputs;
        std::vector<double> inputs;
        for (int i=0; i<memoryOfModel; i++) {
            inputs.push_back(inputsignal[i]);
        }
        
        // 係数alphaの保存
        std::vector<std::vector<double>> A = getPhis(inputsignal, rbfCount, dataCount);
        alpha = solveLeastSquaresMethod(A, inputsignal);
        
        for (int i=0; i<dataCount; i++) {
            double output = 0.0;
            for (int rbfIndex=0; rbfIndex<=rbfCount; rbfIndex++) {
                if (rbfIndex==0) {
                    output += alpha[0];
                }
                else {
                    double squeredNorm = 0.0;
                    for (int j=0; j<memoryOfModel; j++) {
                        squeredNorm += pow((inputsignal[i+j]-rbfs[rbfIndex-1][j]), 2);
                    }
                    output += alpha[rbfIndex] * exp(-spread*squeredNorm);
                }
            }
            
            inputs.erase(inputs.begin());
            inputs.push_back(output);
            
            outputs.push_back(output);
        }
        return outputs;
    }
    
    // k番目の予測関数出力
    std::vector<double> predictionFunctionOutputs(const std::vector<double> &inputsignal, const std::vector<double> &talpha, const int &k, const int &dataCount) const {
        std::vector<double> outputs;
        for (int i=0; i<dataCount; i++) {
            double output = 0.0;
            for (int rbfIndex=0; rbfIndex<=k; rbfIndex++) {
                if (rbfIndex==0) {
                    output += talpha[0];
                }
                else {
                    double squeredNorm = 0.0;
                    for (int j=0; j<memoryOfModel; j++) {
                        squeredNorm += pow((inputsignal[i+j]-rbfs[rbfIndex-1][j]), 2);
                    }
                    output += talpha[rbfIndex] * exp(-spread*squeredNorm);
                }
            }
            outputs.push_back(output);
        }
        return outputs;
    }
    
    // Φを取得する
    std::vector<std::vector<double>> getPhis(const std::vector<double> &inputsignal, const int rbfCount, const int dataCount) const {
        std::vector<std::vector<double>> phis;
        for (int rbfIndex=0; rbfIndex<=rbfCount; rbfIndex++) {
            std::vector<double> vector;
            for (int i=0; i<dataCount; i++) {
                double output;
                if (rbfIndex==0) {
                    output = 1.0;
                }
                else {
                    double squeredNorm = 0.0;
                    for (int j=0; j<memoryOfModel; j++) {
                        squeredNorm += pow((inputsignal[i+j]-rbfs[rbfIndex-1][j]), 2);
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
        
        for (int i=0; i<A.cols(); i++)
            for (int j=0; j<A.rows(); j++)
                A(j, i) = vA[i][j];
        for (int i=0; i<vy.size()-memoryOfModel; i++)
            y(i) = vy[i+memoryOfModel+1];
        
        Eigen::MatrixXd gtg = A.transpose()*A;
        Eigen::MatrixXd gty = A.transpose()*y;
        // Gt*G*o = Gt*y
        Eigen::VectorXd o = (gtg).fullPivHouseholderQr().solve(gty);
        
        std::vector<double> vo(o.size());
        for (int i=0; i<o.size(); i++)
            vo[i] = o[i];
        
        return vo;
    }
    
    int findBestModel(const std::vector<double> &inputsignal) {
        double threshold = 0.1;
        dataCount = inputsignal.size() - memoryOfModel;
        // 最適なモデルが見つかるまでkを増やしながら試行していく
        
        int maxCount = 0;
        int maxRbfIndex = 0;
        for (int rbfIndex=1; rbfIndex<rbfCount; rbfIndex++) {
            std::vector<double> toutput = output(inputsignal, rbfIndex, dataCount);
            std::vector<double> diff;
            for (int i=0; i<dataCount; i++) {
                double diffa = pow(toutput[i] - inputsignal[i], 2);
                diff.push_back(diffa);
            }
            int count = 1;
            while (count < dataCount) {
                double rms = 0;
                for (int i=0; i<count; i++) {
                    rms += diff[i];
                }
                rms /= count;
                rms = sqrt(rms);
                if (rms > threshold) {
                    break;
                }
                
                count++;
            }
            
            if (maxCount < count) {
                maxCount = count;
                maxRbfIndex = rbfIndex;
                std::cout << "[TermRankingMethod] findBestModel - count : " << count << std::endl;
            }
        }
        
        std::cout << "[TermRankingMethod] findBestModel - " << maxRbfIndex << std::endl;
        return maxRbfIndex;
    }
    
    int findBestModelWithError(const std::vector<double> &inputsignal) {
        double threshold = 0.1;
        dataCount = inputsignal.size() - memoryOfModel;
        // 最適なモデルが見つかるまでkを増やしながら試行していく
        
        int maxCount = 0;
        int maxRbfIndex = 0;
        for (int rbfIndex=1; rbfIndex<rbfCount; rbfIndex++) {
            std::vector<double> toutput = output(inputsignal, rbfIndex, dataCount);
            std::vector<double> diff;
            for (int i=0; i<dataCount; i++) {
                double diffa = pow(toutput[i] - inputsignal[i], 2);
                diff.push_back(diffa);
            }
            int count = 1;
            while (count < dataCount) {
                double rms = 0;
                for (int i=0; i<count; i++) {
                    rms += diff[i];
                }
                rms /= count;
                rms = sqrt(rms);
                if (rms > threshold) {
                    break;
                }
                
                count++;
            }
            
            if (maxCount < count) {
                maxCount = count;
                maxRbfIndex = rbfIndex;
                std::cout << "[TermRankingMethod] findBestModel - count : " << count << std::endl;
            }
        }
        
        std::cout << "[TermRankingMethod] findBestModel - " << maxRbfIndex << std::endl;
        return maxRbfIndex;
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
        for (int i=0; i<rbfCount-1; i++) {
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
        for (int i=0; i<rbfCount-1; i++) {
            for (int j=0; j<memoryOfModel; j++) {
                rbfsFstream << rbfs[i][j] << ",";
            }
        }
        rbfsFstream.close();
    }
};

#endif /* defined(__TermRankingMethod__TermRankingMethod__) */
