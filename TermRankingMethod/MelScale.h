//
//  MelScale.h
//  TermRankingMethod
//
//  Created by 狩宿恵介 on 2015/01/29.
//  Copyright (c) 2015年 Keisuke Karijuku. All rights reserved.
//

#ifndef TermRankingMethod_MelScale_h
#define TermRankingMethod_MelScale_h

#include <cmath>


class MelScale {
public:
    
    // Stevens & Volkman 1940; Beranek 1949; O’Shaughnessy 1987
    static double hz2mel_stevens(double f) {
        return 1127.010480 * log(f/700.0 + 1.0);
    }
    
    static double mel2hz_stevens(double mel) {
        return 700.0 * (exp(mel/1127.010480) - 1.0);
    }
    
    
    // Fant 1968
    static double hz2mel_fant(double f) {
        return 1442.695041 * log(f/1000.0 + 1.0);
    }
    
    static double mel2hz_fant(double mel) {
        return 1000.0 * (exp(mel/1442.695041) - 1.0);
    }
    
    
    // Lindsay & Norman 1977
    static double hz2mel_lindsay(double f) {
        return 1046.55994 * log(f/625.0 + 1.0);
    }
    
    static double mel2hz_lindsay(double mel) {
        return 625.0 * (exp(mel/1046.55994) - 1.0);
    }
    
    
    
    /**
     * メルフィルタバンクを作成
     */
    void melFilterBank(double fs, double nfft, int numChannels) {
        
        // ナイキスト周波数
        double fmax = fs / 2;
        
        // ナイキスト周波数
        double melmax = hz2mel_fant(fmax);
        
        // 周波数インデックスの最大値
        double nmax = nfft / 2;
        
        // 周波数解像度（周波数インデックス1あたりのHz幅）
        double df = fs / nfft;
        
        // 各フィルタの中心周波数をお揉める
        double dmel = melmax / (numChannels + 1);
        std::vector<double> melCenters;
        for (int i=1; i<numChannels+1; i++) {
            double melcenter = i * dmel;
            melCenters.push_back(melcenter);
        }
        
        // 各フィルタの中心周波数をHzに変換
        std::vector<double> fCenters;
        for (int i=0; i<numChannels; i++) {
            double hz = mel2hz_fant(melCenters[i]);
            fCenters.push_back(hz);
        }
        
        // 各フィルタの中心周波数を周波数インデックスに変換
        std::vector<double> indexCenters;
        for (int i=0; i<numChannels; i++) {
            double roundedFCenter = round(fCenters[i] / df);
            indexCenters.push_back(roundedFCenter);
        }
        
        // 各フィルタの開始位置のインデックス
        std::vector<double> 
    }
    
    
    
    /*
     
     def melFilterBank(fs, nfft, numChannels):
     """メルフィルタバンクを作成"""
     
     # ナイキスト周波数（Hz）
     fmax = fs / 2
     # ナイキスト周波数（mel）
     melmax = hz2mel(fmax)
     # 周波数インデックスの最大数
     nmax = nfft / 2
     # 周波数解像度（周波数インデックス1あたりのHz幅）
     df = fs / nfft
     # メル尺度における各フィルタの中心周波数を求める
     dmel = melmax / (numChannels + 1)
     melcenters = np.arange(1, numChannels + 1) * dmel
     # 各フィルタの中心周波数をHzに変換
     fcenters = mel2hz(melcenters)
     # 各フィルタの中心周波数を周波数インデックスに変換
     indexcenter = np.round(fcenters / df)
     # 各フィルタの開始位置のインデックス
     indexstart = np.hstack(([0], indexcenter[0:numChannels - 1]))
     # 各フィルタの終了位置のインデックス
     indexstop = np.hstack((indexcenter[1:numChannels], [nmax]))
     
     filterbank = np.zeros((numChannels, nmax))
     for c in np.arange(0, numChannels):
     # 三角フィルタの左の直線の傾きから点を求める
     increment= 1.0 / (indexcenter[c] - indexstart[c])
     for i in np.arange(indexstart[c], indexcenter[c]):
     filterbank[c, i] = (i - indexstart[c]) * increment
     # 三角フィルタの右の直線の傾きから点を求める
     decrement = 1.0 / (indexstop[c] - indexcenter[c])
     for i in np.arange(indexcenter[c], indexstop[c]):
     filterbank[c, i] = 1.0 - ((i - indexcenter[c]) * decrement)
     
     return filterbank, fcenters
     
     
     
     # メルフィルタバンクを作成
     numChannels = 20  # メルフィルタバンクのチャネル数
     df = fs / nfft   # 周波数解像度（周波数インデックス1あたりのHz幅）
     filterbank, fcenters = melFilterBank(fs, nfft, numChannels)
     
     # メルフィルタバンクのプロット
     for c in np.arange(0, numChannels):
     plot(np.arange(0, nfft / 2) * df, filterbank[c])
     savefig("melfilterbank.png")
     show()
     
     */
    
    
    
};

#endif