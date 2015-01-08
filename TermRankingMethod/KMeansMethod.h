//
//  KMeansMethod.h
//  TermRankingMethod
//
//  Created by 狩宿恵介 on 2015/01/08.
//  Copyright (c) 2015年 Keisuke Karijuku. All rights reserved.
//

#ifndef TermRankingMethod_KMeansMethod_h
#define TermRankingMethod_KMeansMethod_h

#include <vector>
#include <random>
#include <limits>

class KMeansMethod {
public:
    
    std::vector<std::vector<double>> calculate(std::vector<double> inputs, int memoryOfModel, double countOfCluster) {
        
        // 最大値、最小値の取得
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::min();
        for (int i=0; i<inputs.size(); i++) {
            if (inputs[i] > max) max = inputs[i];
            if (inputs[i] < min) min = inputs[i];
        }
        
        // クラスタの宣言と初期化
        std::vector<std::vector<double>> cluster;
        std::random_device random_device;
        std::mt19937 mt(random_device());
        std::uniform_real_distribution<double> score(min, max);
        for (int i=0; i<countOfCluster; i++) {
            std::vector<double> vector;
            for (int j=0; j<memoryOfModel; j++) {
                vector.push_back(score(mt));
            }
            cluster.push_back(vector);
        }
        
        // 入力がどこのクラスタに属しているか
        std::vector<int> clusterOfInputs(inputs.size() - memoryOfModel, -1);
        
        while (true) {
            bool isChanged = false;
            
            // 入力をそれぞれ一番近いクラスタに付与
            for (int i=0; i<inputs.size() - memoryOfModel; i++) {
                double minDistance = std::numeric_limits<double>::max();
                double minIndexOfCluster = -1;
                for (int j=0; j<cluster.size(); j++) {
                    double distance = 0;
                    std::vector<double> vector = cluster[j];
                    for (int k=0; k<vector.size(); k++) {
                        distance += pow(inputs[i+k] - vector[k], 2.0);
                    }
                    distance = sqrt(distance);
                    
                    if (distance < minDistance) {
                        minDistance = distance;
                        minIndexOfCluster = j;
                    }
                }
                if (clusterOfInputs[i] != minIndexOfCluster) {
                    isChanged = true;
                    clusterOfInputs[i] = minIndexOfCluster;
                }
            }
            
            if (!isChanged) {
                break;
            }
            
            // 各クラスターから中心を求める
            for (int i=0; i<cluster.size(); i++) {
                std::vector<double> center(memoryOfModel, 0);
                
                std::vector<double> vector = cluster[i];
                std::vector<double> numberOfCluster(countOfCluster, 0);
                for (int j=0; j<memoryOfModel; j++) {
                    for (int k=0; k<inputs.size() - memoryOfModel; k++) {
                        if (clusterOfInputs[k] == i) {
                            center[j] += inputs[k+j];
                            numberOfCluster[j]++;
                        }
                    }
                }
                
                for (int j=0; j<memoryOfModel; j++) {
                    if (numberOfCluster[j] > 0) {
                        center[j] /= numberOfCluster[j];
                    }
                    vector[j] = center[j];
                }
                
            }
        }
        
        return cluster;
    }
};


#endif
