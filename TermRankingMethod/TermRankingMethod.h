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
#include <Eigen/Dense>

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
    void calculate(const std::vector<double> inputsignal) {
        
        KMeansMethod kmeans;
        std::vector<std::vector<double>> cluster = kmeans.calculate(inputsignal, memoryOfModel, rbfCount);
        rbfs = cluster;
        
        std::cout << "[TermRankingMethod] Start Calculate" << std::endl;
        
        dataCount = inputsignal.size() - memoryOfModel;
        
//        // 入力の平均
//        double averageForInput = averageInVector(inputsignal);
//        // 入力の平均と入力の差分
//        double diffInputsignalAndAverage = diffBetweenVectorAndValue(inputsignal, averageForInput);
//        
//        
//        for (int iter=0; iter<rbfCount; iter++) {
//            std::cout << "[TermRankingMethod] Calculate K = " << iter << std::endl;
//            
//            // 線形二乗法による最適化
//            std::vector<std::vector<double>> alphas;
//            for (int k=1; k<rbfCount; k++) {
//                std::vector<std::vector<double>> A = getPhis(inputsignal, k, dataCount);
//                std::vector<double> talpha = solveLeastSquaresMethod(A, inputsignal, k);
//                alphas.push_back(talpha);
//            }
//            
//            
//            // エラー計算
//            std::vector<double> squaredErrors;
//            for (int k=0; k<rbfCount-2; k++) {
//                std::vector<double> predictionOutputs = predictionFunctionOutputs(inputsignal, alphas[k], k, dataCount);
//                double error = errorOfOneStepPrediction(inputsignal, predictionOutputs, diffInputsignalAndAverage);
//                squaredErrors.push_back(error);
//            }
//            
//            
//            // 差分計算
//            std::vector<double> informationDiscrepancys;
//            double firstInformationDiscrepancys = -1.0/inputsignal.size();
//            informationDiscrepancys.push_back(firstInformationDiscrepancys);
//            for (int k=0; k<rbfCount-3; k++) {
//                double informationDiscrepancy = log(squaredErrors[k]/squaredErrors[k+1]);
//                informationDiscrepancys.push_back(informationDiscrepancy);
//            }
//            
//            
//            // 並び替え＆先頭固定
//            for (int tk=iter-1; tk<rbfCount-iter-3; tk++) {
//                for (int k=iter-1; k<rbfCount-iter-3; k++) {
//                    if (informationDiscrepancys[k+1] < informationDiscrepancys[k+2]) {
//                        double tmpA = informationDiscrepancys[k+1];
//                        informationDiscrepancys[k+1] = informationDiscrepancys[k+2];
//                        informationDiscrepancys[k+2] = tmpA;
//                        
//                        auto tmp = rbfs[k+1];
//                        rbfs[k+1] = rbfs[k+2];
//                        rbfs[k+2] = tmp;
//                    }
//                }
//            }
//        }
        
    
        std::cout << "[TermRankingMethod] End Calculate" << std::endl;
    }
    
    std::vector<double> output(const std::vector<double> inputsignal, const int rbfCount, const int dataCount) {
        std::vector<double> outputs;
        std::vector<double> inputs(memoryOfModel);
        for (int i=0; i<memoryOfModel; i++) {
            inputs[i] = inputsignal[i];
        }
        
        // 係数alphaの保存
        std::vector<std::vector<double>> A = getPhis(inputsignal, rbfCount, dataCount);
        alpha = solveLeastSquaresMethod(A, inputsignal, rbfCount);
        
        for (int i=0; i<dataCount; i++) {
            double output = alpha[0];
            for (int rbfIndex=0; rbfIndex<rbfCount-1; rbfIndex++) {
                double squeredNorm = 0.0;
                for (int j=0; j<memoryOfModel; j++) {
                    squeredNorm += pow((inputs[j]-rbfs[rbfIndex][j]), 2);
//                    squeredNorm += pow(inputsignal[i+j]-rbfs[rbfIndex][j], 2);
                }
                output += alpha[rbfIndex] * exp(-spread*squeredNorm);
            }
            
            inputs.erase(inputs.begin());
            inputs.push_back(output);
            
            outputs.push_back(output);
        }
        return outputs;
    }
    
    // k番目の予測関数出力
    std::vector<double> predictionFunctionOutputs(const std::vector<double> inputsignal, const std::vector<double> talpha, const int k, const int dataCount) const {
        std::vector<double> outputs;
        for (int i=0; i<dataCount - memoryOfModel; i++) {
            double output = talpha[0];
            for (int rbfIndex=0; rbfIndex<k-1; rbfIndex++) {
                double squeredNorm = 0.0;
                for (int j=0; j<memoryOfModel; j++) {
                    squeredNorm += pow((inputsignal[i+j]-rbfs[rbfIndex][j]), 2);
                }
                output += talpha[rbfIndex] * exp(-spread*squeredNorm);
            }
            outputs.push_back(output);
        }
        return outputs;
    }
    
    // Φを取得する
    std::vector<std::vector<double>> getPhis(const std::vector<double> inputsignal, const int rbfCount, const int dataCount) const {
        std::vector<std::vector<double>> phis;
        
        std::vector<double> vector;
        for (int j=0; j<rbfCount; j++) {
            vector.push_back(1.0);
        }
        phis.push_back(vector);
        
        for (int i=0; i<dataCount-memoryOfModel-1; i++) {
            std::vector<double> vector;
            for (int rbfIndex=0; rbfIndex<rbfCount; rbfIndex++) {
                double squaredNorm = 0.0;
                for (int k=0; k<memoryOfModel; k++) {
                    squaredNorm += pow(inputsignal[i+k]-rbfs[rbfIndex][k], 2);
                }
                double output = exp(-spread*squaredNorm);

                vector.push_back(output);
            }
            phis.push_back(vector);
        }
        
        return phis;
    }
    
    // ワンステップ誤差を出力
    double errorOfOneStepPrediction(const std::vector<double> inputsignal, const std::vector<double> output, double denominator) const {
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
    
    
    std::vector<double> solveLeastSquaresMethod(const std::vector<std::vector<double>> vA, const std::vector<double> vy, int rbfCount) const {
        Eigen::MatrixXd A((int)(dataCount-memoryOfModel), (int)(rbfCount));
        Eigen::VectorXd y((int)(dataCount-memoryOfModel));
        
        for (int i=0; i<rbfCount; i++)
            for (int j=0; j<dataCount-memoryOfModel; j++)
                A(j, i) = vA[j][i];
        
        for (int i=0; i<dataCount-memoryOfModel; i++)
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
    
    
    int findBestModel(const std::vector<double> inputsignal) {
        double threshold = 0.05;
        dataCount = inputsignal.size() - memoryOfModel;
        // 最適なモデルが見つかるまでkを増やしながら試行していく
        
        int maxCount = 0;
        int maxRbfIndex = 0;
        for (int rbfIndex=1; rbfIndex<rbfCount; rbfIndex++) {
            std::vector<double> toutput = output(inputsignal, rbfIndex, dataCount);
            std::vector<double> diff;
            for (int i=0; i<dataCount - memoryOfModel; i++) {
                double tmp = pow(toutput[i] - inputsignal[i+1], 2);
                diff.push_back(tmp);
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
    
    int findBestModelWithError(const std::vector<double> inputsignal) {
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
    
    double averageInVector(const std::vector<double> vector) const {
        double average = 0.0;
        auto iter = vector.begin();
        auto iter_end = vector.end();
        while (iter != iter_end) {
            average += (*iter);
            ++iter;
        }
        return average/vector.size();
    }
    
    double diffBetweenVectorAndValue(const std::vector<double> vector, const double value) const {
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
