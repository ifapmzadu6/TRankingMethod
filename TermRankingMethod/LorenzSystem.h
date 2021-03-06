//
//  LorenzSystem.h
//  TermRankingMethod
//
//  Created by Keisuke Karijuku on 2014/12/20.
//  Copyright (c) 2014年 Keisuke Karijuku. All rights reserved.
//

#ifndef TermRankingMethod_LorenzSystem_h
#define TermRankingMethod_LorenzSystem_h

class LorenzSystem {
public:
    
    double x;
    double y;
    double z;
    
    double p = 10.0;
    double r = 28.0;
    double b = 3.0/8.0;
    
    double dt = 0.02;
    
    LorenzSystem() {
        x = 3.584036031;
        y = 1.398178981;
        z = 25.08707469;
//        x = 1;
//        y = 1;
//        z = 1;
    }
    
    void nextTime() {
        x += dx(x, y, z);
        y += dy(x, y, z);
        z += dz(x, y, z);
    }
    
    void nextTimeWithRungetta() {
        double k1_dx = dx(x, y, z);
        double k1_dy = dy(x, y, z);
        double k1_dz = dz(x, y, z);
        
        double k2_dx = dx( x+k1_dx*0.5, y+k1_dy*0.5, z+k1_dz*0.5 );
        double k2_dy = dy( x+k1_dx*0.5, y+k1_dy*0.5, z+k1_dz*0.5 );
        double k2_dz = dz( x+k1_dx*0.5, y+k1_dy*0.5, z+k1_dz*0.5 );
        
        double k3_dx = dx( x+k2_dx*0.5, y+k2_dy*0.5, z+k2_dz*0.5 );
        double k3_dy = dy( x+k2_dx*0.5, y+k2_dy*0.5, z+k2_dz*0.5 );
        double k3_dz = dz( x+k2_dx*0.5, y+k2_dy*0.5, z+k2_dz*0.5 );
        
        double k4_dx = dx( x+k3_dx, y+k3_dy, z+k3_dz );
        double k4_dy = dy( x+k3_dx, y+k3_dy, z+k3_dz );
        double k4_dz = dz( x+k3_dx, y+k3_dy, z+k3_dz );
        
        x += ( k1_dx + 2.0*k2_dx + 2.0*k3_dx + k4_dx ) / 6.0;
        y += ( k1_dy + 2.0*k2_dy + 2.0*k3_dy + k4_dy ) / 6.0;
        z += ( k1_dz + 2.0*k2_dz + 2.0*k3_dz + k4_dz ) / 6.0;
    }
    
    double dx(double x, double y, double z) {
        double dx = (-p * x + p * y) * dt;
        return dx;
    }
    
    double dy(double x, double y, double z) {
        double dy = (-x * z + r * x - y) * dt;
        return dy;
    }
    
    double dz(double x, double y, double z) {
        double dz = (x * y - b * z) * dt;
        return dz;
    }
    
};

#endif
