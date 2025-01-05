#pragma once

#include <memory>

#include "nlohmann/json.hpp"
#include "teqp/cpp/teqpcpp.hpp"

namespace teqpflsh::properties::interfaces{

using namespace teqp::cppinterface;

// ABC for interfaces, defining what functions must be implemented by thermodynamic model implementations
class HelmholtzInterface{
public:
    virtual Eigen::Array33d get_Armat2(const double T, const double rho, const Eigen::ArrayXd& z) const = 0;
    virtual std::tuple<double, double, double> get_A00A10A01(const double T, const double rho, const Eigen::ArrayXd& z) const = 0;
    virtual ~HelmholtzInterface() = default;
};

class teqpHelmholtzInterface : public HelmholtzInterface{
public:
    std::unique_ptr<AbstractModel> model_ideal_gas, model_residual;
    
    teqpHelmholtzInterface(const std::string &ideal_gas, const std::string &resid)
        : model_ideal_gas(make_model(nlohmann::json::parse(ideal_gas))),
          model_residual(make_model(nlohmann::json::parse(resid))) {}
    
    std::tuple<double, double, double> get_A00A10A01(const double T, const double rho, const Eigen::ArrayXd& z) const override {
        Eigen::Array2d A00A01 = model_ideal_gas->get_Ar01n(T, rho, z) + model_residual->get_Ar01n(T, rho, z); // This gives [A00, A01] at the same time
        auto A10 = model_ideal_gas->get_Ar10(T, rho, z) + model_residual->get_Ar10(T, rho, z);
        return std::make_tuple(A00A01[0], A10, A00A01[1]);
    };
    
    Eigen::Array33d get_Armat2(double T, double rho, const Eigen::ArrayXd& z) const override{
        return model_ideal_gas->get_deriv_mat2(T, rho, z) + model_residual->get_deriv_mat2(T, rho, z);
    }
};

}
