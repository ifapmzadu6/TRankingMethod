//
//  main.cpp
//  TermRankingMethod
//
//  Created by Keisuke Karijuku on 2014/12/01.
//  Copyright (c) 2014年 Keisuke Karijuku. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "TermRankingMethod.h"

#include "Wave.h"
//#include "analysis.h"


int main(int argc, const char * argv[]) {
    
    int rbfCount = 150;
    int memoryOfModel = 8;
    
    int dataCount = 2000;
    int startDataCount = 1000;
    
    // Get sound wave.
    std::vector<double> tmp;
    Wave wav;
    if(wav.InputWave("/Users/KeisukeKarijuku/Dropbox/sample.wav") != 0)
        return -1;
    wav.StereoToMono();
    wav.GetData(tmp);
    std::cout << tmp.size() << std::endl;
    
    // Nomilization
    double max = 0, min = std::numeric_limits<double>::max();
    for (int i=startDataCount; i<startDataCount+dataCount; i++) {
        if (max < tmp[i]) max = tmp[i];
        if (min > tmp[i]) min = tmp[i];
    }
    for (int i=startDataCount; i<startDataCount+dataCount; i++) {
        tmp[i] = (tmp[i] - min)/(max - min);
    }
    
    // Cropped Input
    std::vector<double> inputsignal;
    for (int i=startDataCount; i<startDataCount+dataCount; i++) {
        inputsignal.push_back(tmp[i]);
    }
    
    // Term Ranking Method
    TermRankingMethod termRankingMethod(rbfCount, memoryOfModel);
    termRankingMethod.calculate(inputsignal);
    std::vector<double> outputsignal = termRankingMethod.output(inputsignal, dataCount);
    
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
