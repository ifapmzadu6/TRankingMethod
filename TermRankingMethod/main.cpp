//
//  main.cpp
//  TermRankingMethod
//
//  Created by Keisuke Karijuku on 2014/12/01.
//  Copyright (c) 2014年 Keisuke Karijuku. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <limits>

#include "TermRankingMethod.h"

#include "Wave.h"
//#include "analysis.h"

#include "LorenzSystem.h"



int main(int argc, const char * argv[]) {
    
    int rbfCount = 400;
    int memoryOfModel = 7;
    
    int dataCount = 3000;
    int startDataCount = 1000;
    
    bool readFromFile = false;
    
    // Get sound wave.
    std::vector<double> tmp;
    Wave wav;
    if(wav.InputWave("sample.wav") != 0)
        return -1;
    wav.StereoToMono();
    wav.Normalize();
    wav.GetData(tmp);
    std::cout << tmp.size() << std::endl;
    
    // Cropped Input
    std::vector<double> inputsignal;
    
    for (int i=startDataCount; i<startDataCount+dataCount; i++) {
        double signal = (tmp[i] + 1.0) / 2.0;
        inputsignal.push_back(signal);
    }
    
    std::ofstream inputStream("inputs.txt");
    for (int i=0; i<inputsignal.size(); i++) {
        inputStream << inputsignal[i] << std::endl;
    }
    inputStream.close();
    
//    LorenzSystem lorenzSystem = LorenzSystem();
//    std::ofstream lorenzSystemStream("lorenzSystem.txt");
//    for (int i=0; i<dataCount; i++) {
//        lorenzSystemStream << lorenzSystem.x << " " << lorenzSystem.y << " " << lorenzSystem.z << std::endl;
//        
//        inputsignal.push_back(lorenzSystem.x);
//        
//        lorenzSystem.nextTime();
//    }
//    
//    double max = std::numeric_limits<double>::min();
//    double min = std::numeric_limits<double>::max();
//    for (int i=0; i<inputsignal.size(); i++) {
//        if (inputsignal[i] > max) {
//            max = inputsignal[i];
//        }
//        if (inputsignal[i] < min) {
//            min = inputsignal[i];
//        }
//    }
    
    // Term Ranking Method
    TermRankingMethod termRankingMethod(rbfCount, memoryOfModel);
    if (readFromFile) {
        termRankingMethod.readAlpha("alpha.txt");
        termRankingMethod.readRBFs("rbfs.txt");
    }
    else {
        termRankingMethod.calculate(inputsignal);
    }
    
    int bestRbfCount = termRankingMethod.findBestModel(inputsignal);
    std::vector<double> outputsignal = termRankingMethod.output(inputsignal, rbfCount, dataCount - memoryOfModel);
    
    termRankingMethod.writeAlpha("alpha.txt");
    termRankingMethod.writeRBFs("rbfs.txt");
    
    
    // Error
    double error = 0;
    for (int i=0; i<dataCount; i++) {
        error += pow((inputsignal[i]-outputsignal[i]), 2);
    }
    std::cout << error << std::endl;
    
    
    // Gnuplot
    std::ofstream fstream("result.txt");
    for (int i=0; i<dataCount; i++) {
        fstream << i << " ";
        fstream << inputsignal[i] << " ";
        fstream << outputsignal[i] << std::endl;
    }
    fstream.close();
    std::system("/usr/local/bin/gnuplot -persist -e \" p 'result.txt' u 1:2 w l, '' u 1:3 w l \"");
    
    return 0;
}
