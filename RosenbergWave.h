//
//  RosenbergWave.h
//  TermRankingMethod
//
//  Created by 狩宿恵介 on 2015/02/01.
//  Copyright (c) 2015年 Keisuke Karijuku. All rights reserved.
//

#ifndef TermRankingMethod_RosenbergWave_h
#define TermRankingMethod_RosenbergWave_h

class RosenbergWave {
public:
    
    /**
     * Rosenberg波を生成
     */
    static double GenRosenberg(int freq, double samplePerSec) {
        static double t;
        double tau = 0.90;    /* 声門開大期 */
        double tau2 = 0.95;   /* 声門閉小期 */
        double sample = 0.0;
        
        t += (double)freq / samplePerSec;
        t -= floor(t);
        
        if (t <= tau) {
            sample = 3.0*pow(t/tau,2.0)-2.0*pow(t/tau,3.0);
        } else if (t < tau+tau2) {
            sample = 1.0-pow((t-tau)/tau2,2.0);
        }
        
        return 2.0*(sample-0.5);
    }
};

#endif
