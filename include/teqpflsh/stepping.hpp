#pragma once

#include <memory>
#include <Eigen/Dense>

#include <iostream>

namespace teqpflsh{

template<typename M>
auto build_rJ(const M& ALPHA, double T, double rho, double R, char x_key, double x_val, char y_key, double y_val)
{

    Eigen::Vector2d r; 
    Eigen::Matrix2d J;
    double w = 1/T, dwdT = -w*w;

    for (auto i = 0; i < 2; ++i){
        char key = (i == 0) ? x_key : y_key;
        double target = (i == 0) ? x_val : y_val;
        switch(key){
            case 'S':{
                auto w2dadw = R*(ALPHA(1,0) - ALPHA(0,0));
                auto wdadw = w2dadw/w;

                auto w2d2adw2 = R/w*(ALPHA(2,0) - 2*ALPHA(1,0) + 2*ALPHA(0,0));
                auto w2d2adwdrho = R*(ALPHA(1,1)/rho - ALPHA(0,1)/rho);
                
                r(i) = w2dadw - target;

                // ds/dT|rho = ds/d(1/T)|rho*d(1/T)/dT
                J(i,0) = (w2d2adw2 + 2*wdadw)*dwdT;
                J(i,1) = w2d2adwdrho; // ds/drho|T
                break;
            }
            case 'H':{
                auto w2dadw = R*(ALPHA(1,0) - ALPHA(0,0));
                auto wdadw = w2dadw/w;
                auto dadw = wdadw/w;

                auto rhodadrho = R/w*ALPHA(0,1);
                auto dadrho = rhodadrho/rho;
                auto rhod2adwdrho = R/(w*w)*(ALPHA(1,1) - ALPHA(0,1));
                auto rhod2adrho2 = R/w*ALPHA(0,2)/rho;
                auto w2d2adw2 = R/w*(ALPHA(2,0) - 2*ALPHA(1,0) + 2*ALPHA(0,0));
                auto wd2adw2 = w2d2adw2/w;
                auto w2d2adwdrho = R*(ALPHA(1,1)/rho - ALPHA(0,1)/rho);
                auto wd2adwdrho = w2d2adwdrho/w;

                r(i) = (ALPHA(0,0)*R*T + wdadw + rhodadrho) - target;
                // dh/dT|rho = dh/d(1/T)|rho*d(1/T)/dT
                J(i,0) = (wd2adw2 + rhod2adwdrho + 2*dadw)*dwdT;
                J(i,1) = wd2adwdrho + rhod2adrho2 + 2*dadrho; // dh/drho|T
                break;
            }
            default:
                throw std::invalid_argument("");
        }
    }
    return std::make_tuple(r, J);
}

// class Newton2DStepper{
//     public:
//         std::unique_ptr<teqp::AbstractModel> model_ig, model_res;
//         char x_key, y_key;

//         Newton2DStepper(model_ig, model_res, x_key, x_val, y_key, y_val) {};

//         auto get_rJ(x){
//             auto [T, rho] = x;
//             auto dh_ig = model_ig->DerivativeHolderSquare<2>(T, rho, z);
//             auto dh_res = model_res->DerivativeHolderSquare<2>(T, rho, z);
//             auto ALPHA = dh_ig + dh_res;
//             auto [r,J] = build_rJ(ALPHA);
//         }
// };

} /* namespace teqpflsh */