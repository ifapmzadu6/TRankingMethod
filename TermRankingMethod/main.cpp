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
#include "analysis.h"
#include "LorenzSystem.h"
#include "Audio.h"
#include "Window.h"
#include "MelScale.h"


int main(int argc, const char * argv[]) {
    
    using namespace std;
    
    MelScale a;
    vector<FilterBank> melFilterBank = a.melFilterBank(44100, 44100, 20);
    
    MelScale::plotMelFilterBank(44100, 44100, 20);
    
    int rbfCount = 100;
    int memoryOfModel = 7;
    double spread = 1.0 / (2*1.);
    int dataCount = 100000;
    int startDataCount = 5000;
    
//    int rbfCount = 150;
//    int memoryOfModel = 2;
//    int dataCount = 1000;
//    int startDataCount = 0;
//    double spread = 1.0 / (2*1.);
    
    bool readFromFile = false;
    
    // Get sound wave.
    vector<double> tmp;
    vector<int> fp;
    Wave wav;
    if(wav.InputWave("sample.wav") != 0)
        return -1;
    wav.StereoToMono();
    wav.Normalize();
    wav.GetData(tmp);
    GetFlucPeriod(fp, tmp);
    cout << fp.size() << endl;
    
    vector<double> inputsignal;
    
    
    for (int i=startDataCount; i<startDataCount+dataCount; i++) {
        double signal = (tmp[i] + 1.0) / 2.0;
        inputsignal.push_back(signal);
    }
    
    vector<double> windowed = Window::hamming(inputsignal);
    
    vector<AudioComplex> outputs = Audio::dft_r2c_1d_vector(windowed, 44100, 0);
    
    
    
//     Nomilization
//    double max = 0, min = numeric_limits<double>::max();
//    for (int i = startDataCount; i < startDataCount + dataCount; i++) {
//        if (max < fp[i]) max = fp[i];
//        if (min > fp[i]) min = fp[i];
//    }
//    for (int i = startDataCount; i < startDataCount + dataCount; i++) {
//        double value = ((double)fp[i] - min) / (max - min);
//        inputsignal.push_back(value);
//    }
    
//    // ただのテストデータ
//    for (int i=0; i<dataCount; i++) {
//        double value = 1.0 * i/dataCount;
//        inputsignal.push_back(value);
//    }
    
    
//    // ローレンツシステム
//    LorenzSystem lorenzSystem = LorenzSystem();
//    ofstream lorenzSystemStream("lorenzSystem.txt");
//    for (int i=0; i<dataCount; i++) {
//        lorenzSystemStream << lorenzSystem.x << " " << lorenzSystem.y << " " << lorenzSystem.z << endl;
//        
//        inputsignal.push_back(lorenzSystem.x);
//        
//        lorenzSystem.nextTime();
//    }
//    
//    double max = numeric_limits<double>::min();
//    double min = numeric_limits<double>::max();
//    for (int i=0; i<inputsignal.size(); i++) {
//        if (inputsignal[i] > max) {
//            max = inputsignal[i];
//        }
//        if (inputsignal[i] < min) {
//            min = inputsignal[i];
//        }
//    }
    
    vector<double> amp_dft;
    for (int i=0; i<outputs.size(); i++) {
        double amp = sqrt(pow(outputs[i].re, 2) + pow(outputs[i].im, 2));
        double log_amp = log10(amp);
        amp_dft.push_back(log_amp);
    }
    
    double index = 0.0;
    ofstream fft("fft.txt");
    while (true) {
        double mel = MelScale::mel2hz_stevens(index);
        if (index >= outputs.size()) {
            break;
        }
        fft << mel << " ";
        fft << amp_dft[index] << endl;
        
        index++;
    }
    fft.close();
    
    system("/usr/local/bin/gnuplot -persist -e \"set xr [0:5000]; p 'fft.txt' u 1:2 w l \"");
    

//    // Term Ranking Method
//    TermRankingMethod termRankingMethod(rbfCount, memoryOfModel, spread);
//    if (readFromFile) {
//        termRankingMethod.readAlpha("alpha.txt");
//        termRankingMethod.readRBFs("rbfs.txt");
//    }
//    else {
//        termRankingMethod.calculate(inputsignal);
//    }
//    
////    int bestRbfCount = termRankingMethod.findBestModel(inputsignal);
////    vector<double> outputsignal = termRankingMethod.output(inputsignal, bestRbfCount, dataCount);
//    vector<double> outputsignal = termRankingMethod.output(inputsignal, rbfCount, dataCount);
//    
//    termRankingMethod.writeAlpha("alpha.txt");
//    termRankingMethod.writeRBFs("rbfs.txt");
//    
//    
//    // Error
//    double error = 0;
//    for (int i=0; i<dataCount; i++) {
//        error += pow((inputsignal[i]-outputsignal[i]), 2);
//    }
//    cout << error << endl;
//    
//    
//    // Gnuplot
//    ofstream fstream("result.txt");
//    for (int i=0; i<dataCount; i++) {
//        fstream << i << " ";
//        fstream << inputsignal[i] << " ";
//        fstream << outputsignal[i] << endl;
//    }
//    fstream.close();
//    system("/usr/local/bin/gnuplot -persist -e \"set xr [0:1000]; p 'result.txt' u 1:2 w l, '' u 1:3 w l \"");
//    
//    // 環境を保存
//    ofstream tfstream("環境.txt");
//    tfstream << "rbfCount = " << rbfCount << endl;
//    tfstream << "memoryOfModel = " << memoryOfModel << endl;
//    tfstream << "spread = " << spread << endl;
//    tfstream << "dataCount = " << dataCount << endl;
//    tfstream << "startDataCount = " << startDataCount << endl;
//    tfstream.close();
    
    
    return 0;
}









