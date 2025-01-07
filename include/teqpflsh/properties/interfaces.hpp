#pragma once

#include <memory>

#include "nlohmann/json.hpp"
#include "teqp/cpp/teqpcpp.hpp"

namespace teqpflsh::properties::interfaces{

using namespace teqp::cppinterface;

/** ABC for interfaces, defining what functions must be implemented by thermodynamic model implementations

 We use the nomenclature
\f[
\Lambda_{ij} = \varpi^i\rho^j\deriv{^{i+j}(\alpha)}{\varpi^i\partial\rho^j}{}
\f]
 with \f$\varpi=1/T\f$
*/
class HelmholtzInterface{
public:
    /**
     Return the matrix of values obtained from derivatives of the reduced Helmholtz energy. Entries are of the form
     \f[
     A^{\rm r}_{i,j} =\Lambda^{\rm (ig)}_{ij} + \Lambda^{\rm r}_{ij}
     \f]
     where it is the sum of the residual and ideal-gas portions that go into the matrix
     
     The matrix is constructed in this way because autodiff allows the rows and columns of pure derivatives (i=0 and j=0) to be obtained as a single operation
     \param T Temperature, in K
     \param rho Molar density, in mol/m^3
     \param molefrac mole fractions
     */
    virtual Eigen::Array33d get_Armat2(const double T, const double rho, const Eigen::ArrayXd& molefrac) const = 0;
    /**
     \brief Return a tuple of the values \f$\Lambda_{00}\f$, \f$\Lambda_{10}\f$, and \f$\Lambda_{01}\f$ where each term is a sum of the residual and ideal-gas contributions
     \param T Temperature, in K
     \param rho Molar density, in mol/m^3
     \param molefrac mole fractions
     */
    virtual std::tuple<double, double, double> get_A00A10A01(const double T, const double rho, const Eigen::ArrayXd& molefrac) const = 0;
    
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
