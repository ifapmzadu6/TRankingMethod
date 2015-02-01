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
#include "RosenbergWave.h"


std::vector<double> getMFCC(std::vector<double> inputsignal, double samplePerSec);




int main(int argc, const char * argv[]) {
    
    using namespace std;
    
    int rbfCount = 100;
    int memoryOfModel = 1;
    double spread = 1.0 / (2*0.0001);
    int dataCount = 100000;
    int startDataCount = 5000;
    
//    int rbfCount = 150;
//    int memoryOfModel = 2;
//    int dataCount = 1000;
//    int startDataCount = 0;
//    double spread = 1.0 / (2*1.);
    
    bool readFromFile = false;
    
    // Get sound wave.
    vector<double> wave;
    vector<int> fp;
    Wave wav;
    if(wav.InputWave("sample.wav") != 0)
        return -1;
    wav.StereoToMono();
    wav.Normalize();
    wav.GetData(wave);
    int samplePerSec = wav.GetSamplesPerSec();
    
    // Analyze
    GetFlucPeriod(fp, wave);
    cout << fp.size() << endl;
    
    vector<double> inputsignal;
    
    
    for (int i=startDataCount; i<startDataCount+dataCount; i++) {
        double signal = (wave[i] + 1.0) / 2.0;
        inputsignal.push_back(signal);
    }
    
    vector<double> windowed = Window::hamming(inputsignal);
    
    double fft_size = 44100;
    vector<AudioComplex> outputs = Audio::dft_r2c_1d_vector(windowed, fft_size, 0);
    
    vector<double> voice_mfcc = getMFCC(inputsignal, samplePerSec);
    
    vector<double> rosenberg_input;
    for (int i=0; i<dataCount; i++) {
        double signal = (RosenbergWave::GenRosenberg(20000, 44100) + 1.0) / 2.0;
        rosenberg_input.push_back(signal);
    }
    
    vector<double> rosenberg_mfcc = getMFCC(rosenberg_input, samplePerSec);
    
    
    
    
    
    
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
//
    
    
    
    
    
    
    
    
    
    
    // Term Ranking Method
    TermRankingMethod termRankingMethod(rbfCount, memoryOfModel, spread);
    if (readFromFile) {
        termRankingMethod.readAlpha("alpha.txt");
        termRankingMethod.readRBFs("rbfs.txt");
    }
    else {
        termRankingMethod.calculate(voice_mfcc);
    }
    
//    int bestRbfCount = termRankingMethod.findBestModel(voice_mfcc);
    //    vector<double> outputsignal = termRankingMethod.output(inputsignal, bestRbfCount, dataCount);
    vector<double> outputsignal = termRankingMethod.output(voice_mfcc, rbfCount, (int)voice_mfcc.size());
    
    termRankingMethod.writeAlpha("alpha.txt");
    termRankingMethod.writeRBFs("rbfs.txt");
    
    
    // Error
    double error = 0;
    for (int i=0; i<voice_mfcc.size(); i++) {
        error += pow((voice_mfcc[i]-outputsignal[i]), 2);
    }
    cout << error << endl;
    
    
    // Gnuplot
    ofstream fstream("result.txt");
    for (int i=0; i<voice_mfcc.size(); i++) {
        fstream << i << " ";
        fstream << voice_mfcc[i] << " ";
        fstream << outputsignal[i] << endl;
    }
    fstream.close();
    system("/usr/local/bin/gnuplot -persist -e \"p 'result.txt' u 1:2 w l, '' u 1:3 w l \"");
    
    // 環境を保存
    ofstream tfstream("環境.txt");
    tfstream << "rbfCount = " << rbfCount << endl;
    tfstream << "memoryOfModel = " << memoryOfModel << endl;
    tfstream << "spread = " << spread << endl;
    tfstream << "dataCount = " << dataCount << endl;
    tfstream << "startDataCount = " << startDataCount << endl;
    tfstream.close();
    
    
    
    
    
    return 0;
}









std::vector<double> getMFCC(std::vector<double> inputsignal, double samplePerSec) {
    
    using namespace std;
    
    vector<double> windowed = Window::hamming(inputsignal);
    
    double fft_size = 44100;
    vector<AudioComplex> outputs = Audio::dft_r2c_1d_vector(windowed, fft_size, 0);
    
    
    vector<double> amp_dft;
    for (int i=0; i<fft_size/2; i++) {
        double amp = sqrt(pow(outputs[i].re, 2) + pow(outputs[i].im, 2));
        amp_dft.push_back(amp);
    }
    
    double index = 0.0;
    
    ofstream fft("fft.txt");
    while (true) {
        double mel = MelScale::mel2hz_stevens(index);
        if (mel >= amp_dft.size()) {
            break;
        }
        fft << mel << " ";
        fft << log10(amp_dft[mel]) << endl;
        
        index++;
    }
    fft.close();
    
//    system("/usr/local/bin/gnuplot -persist -e \"p 'fft.txt' w l \"");
    
    // メルフィルタバンクに変換
    ofstream melfilterbank("melfilterbank.txt");
    
    vector<MelFilterBank> melFilterBank = MelScale::melFilterBank(samplePerSec, fft_size, 96);
    vector<double> melSpectrum;
    for (int i=0; i<melFilterBank.size(); i++) {
        double d = 0;
        for (int j=0; j<(melFilterBank[i].stopIndex - melFilterBank[i].startIndex); j++) {
            double s = melFilterBank[i].filter[j] * amp_dft[j + melFilterBank[i].startIndex];
            d += s;
        }
        melfilterbank << melFilterBank[i].centerIndex << " ";
        melfilterbank << log10(d) << endl;
        
        melSpectrum.push_back(log10(d));
    }
    
    melfilterbank.close();
    
//    system("/usr/local/bin/gnuplot -persist -e \"p 'melfilterbank.txt' w l\"");
    
    
    // メル周波数ケプストラム
    vector<double> dct = Audio::dct_r2r_1d_vector(melSpectrum);
    
    ofstream melCepstral("melcepstral.txt");
    
    for (int i=0; i<dct.size(); i++) {
        melCepstral << dct[i] << endl;
    }
    
    melCepstral.close();
    
//    system("/usr/local/bin/gnuplot -persist -e \"p 'melcepstral.txt' w l\"");
    
    
    // メル周波数ケプストラム係数(0次元を除き、mfcc_N次元だけ利用する)
    int mfcc_N = 12;
    vector<double> mfcc;
    for (int i=1; i<mfcc_N+1; i++) {
        mfcc.push_back(dct[i]);
    }
    
    // MFCCの平均除去
    double ave = 0;
    for (int i=0; i<mfcc.size(); i++) {
        ave += mfcc[i];
    }
    ave /= mfcc.size();
    for (int i=0; i<mfcc.size(); i++) {
        mfcc[i] = mfcc[i] - ave;
    }
    
    // MFCCのノーマライズ
    double abs_max = 0;
    for (int i=0; i<mfcc.size(); i++) {
        if (abs(mfcc[i]) > abs_max) {
            abs_max = mfcc[i];
        }
    }
    for (int i=0; i<mfcc.size(); i++) {
        mfcc[i] = mfcc[i] / abs_max;
    }
    
    
    ofstream mfccfstream("mfcc.txt");
    
    for (int i=0; i<mfcc.size(); i++) {
        mfccfstream << mfcc[i] << endl;
    }
    
    mfccfstream.close();
    
    system("/usr/local/bin/gnuplot -persist -e \"p 'mfcc.txt' w l\"");
    
    return mfcc;
}
