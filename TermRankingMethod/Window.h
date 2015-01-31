//
//  Window.h
//  TermRankingMethod
//
//  Created by 狩宿恵介 on 2015/01/26.
//  Copyright (c) 2015年 Keisuke Karijuku. All rights reserved.
//

#ifndef TermRankingMethod_Window_h
#define TermRankingMethod_Window_h

#include <vector>


namespace Window {
    
    std::vector<double> hamming(std::vector<double> input) {
        std::vector<double> output;
        for (int i=0; i<input.size(); i++) {
            double window = 0.54 - 0.46 * cos(2 * M_PI * input[i]);
            output.push_back(window);
        }
        return output;
    }
    
}

#endif
