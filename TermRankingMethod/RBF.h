//
//  RBF.h
//  FireflyProject
//
//  Created by Keisuke Karijuku on 2014/09/22.
//  Copyright (c) 2014年 Keisuke Karijuku. All rights reserved.
//

#ifndef __FireflyProject__RBF__
#define __FireflyProject__RBF__

#include <cmath>
#include <vector>

double function(const double &spread, const std::vector<double> &centerVector, const std::vector<double> &x);
double mse(const std::vector<std::vector<double>> &d, const std::vector<std::vector<double>> &o);
double squeredNorm(const std::vector<double> &a, const std::vector<double> &b);
void mult(std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);

#endif /* defined(__FireflyProject__RBF__) */
