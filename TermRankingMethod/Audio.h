//
//  Audio.h
//  TermRankingMethod
//
//  Created by 狩宿恵介 on 2015/01/25.
//  Copyright (c) 2015年 Keisuke Karijuku. All rights reserved.
//

#ifndef TermRankingMethod_Audio_h
#define TermRankingMethod_Audio_h

#include <iostream>
#include <fstream>
#include <vector>
#include <fftw3.h>


struct AudioComplex {
    double re;
    double im;
};


class Audio {
public:
    
    
    
    /**
     * 実データ1次元離散フーリエ変換（One-Dimensioanal DFTs of Real Data）
     */
    static std::vector<AudioComplex> dft_r2c_1d_vector(std::vector<double> input_vector, int fft_size, unsigned flags) {
        int n = (int)input_vector.size();
        double *input = (double *)malloc(sizeof(double) * n);
        for (int i=0; i<n; i++) {
            input[i] = input_vector[i];
        }
        fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
        
        fftw_plan plan = fftw_plan_dft_r2c_1d(fft_size, input, output, FFTW_ESTIMATE);
        fftw_execute(plan);
        
        std::vector<AudioComplex> output_vector;
        for (int i=0; i<n/2+1; i++) {
            AudioComplex complex = {output[i][0], output[i][1]};
            output_vector.push_back(complex);
        }
        
        if(plan) fftw_destroy_plan(plan);
        fftw_free(output);
        free(input);
        
        return output_vector;
    }
    
    
    
    /**
     * 実データ1次元離散コサイン変換
     */
    static std::vector<double> dct_r2r_1d_vector(std::vector<double> input_vector) {
        int n = (int)input_vector.size();
        double *input = (double *)malloc(sizeof(double) * n);
        double *output = (double *)malloc(sizeof(double) * n);
        
        for (int i=0; i<n; i++) {
            input[i] = input_vector[i];
        }
        
        fftw_plan plan = fftw_plan_r2r_1d(n, input, output, FFTW_REDFT10, FFTW_ESTIMATE);
        fftw_execute(plan);
        
        std::vector<double> output_vector;
        for (int i=0; i<n; i++) {
            output_vector.push_back(output[i]);
        }
        
        if(plan) fftw_destroy_plan(plan);
        free(output);
        free(input);
        
        return output_vector;
    }
    
    /**
     * 実データ1次元離散コサイン逆変換
     */
    static std::vector<double> dct_r2r_1d_reverse_vector(std::vector<double> input_vector) {
        int n = (int)input_vector.size();
        double *input = (double *)malloc(sizeof(double) * n);
        double *output = (double *)malloc(sizeof(double) * n);
        
        for (int i=0; i<n; i++) {
            input[i] = input_vector[i];
        }
        
        fftw_plan plan = fftw_plan_r2r_1d(n, input, output, FFTW_REDFT01, FFTW_ESTIMATE);
        fftw_execute(plan);
        
        std::vector<double> output_vector;
        for (int i=0; i<n; i++) {
            output_vector.push_back(output[i]);
        }
        
        if(plan) fftw_destroy_plan(plan);
        free(output);
        free(input);
        
        return output_vector;
    }
};


#endif
