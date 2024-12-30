#define NOMINMAX

#include <cstdlib>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
using Catch::Matchers::WithinAbsMatcher;
using Catch::Matchers::WithinRelMatcher;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

#include <numeric>
#include <map>
#include <string>
#include <set>
#include <filesystem>

#include "teqp/json_tools.hpp"
#include "teqp/ideal_eosterms.hpp"
#include "teqpflsh/superancillary/superancillary.hpp"
#include "teqpflsh/properties/interfaces.hpp"

#include<mach/mach.h>
auto get_resident_memory(){
    
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    
    if (KERN_SUCCESS != task_info(mach_task_self(),
                                  TASK_BASIC_INFO, (task_info_t)&t_info,
                                  &t_info_count))
    {
        return -1UL;
    }
    else{
        return t_info.resident_size;
    }
};

// A class storing the info for a single two-phase point
struct TwoPhasePoint{
    double T;
    double Q;
    double D;
    double H;
    double S;
    double U;
    double P;
    
    auto get_val(char k){
        switch (k){
            case 'T': return T;
            case 'D': return D;
            case 'P': return P;
            case 'Q': return Q;
            case 'H': return H;
            case 'S': return S;
            case 'U': return U;
            default: throw std::invalid_argument(std::string(1,k));
        }
    }
};

struct TwoPhaseResult{
    double Terr;
    double Qerr;
    double elap_us;
    double count;
//    stdproppair;
};

using namespace teqpflsh::superancillary;
using namespace teqpflsh::properties;

class SuperancillaryFixture{
private:
    auto build_points(auto Tvec, auto qvec){
        std::vector<TwoPhasePoint> points;
        for (auto T : Tvec){
            for (auto q : qvec){
                TwoPhasePoint pt;
                pt.T = T;
                pt.Q = q;
                pt.P = sa.get_yval(T, q, 'P');
                pt.D = sa.get_yval(T, q, 'D');
                pt.H = sa.get_yval(T, q, 'H');
                pt.S = sa.get_yval(T, q, 'S');
                pt.U = sa.get_yval(T, q, 'U');
                points.push_back(pt);
            }
        }
        return points;
    }
public:
    SuperAncillary<Eigen::ArrayXd> sa;
    std::vector<TwoPhasePoint> points;
    
    SuperancillaryFixture() : sa(SuperAncillary(teqp::load_a_JSON_file("../src/testdata/WATER_exps.json")))
    {
        nlohmann::json FLD = teqp::load_a_JSON_file("../externals/teqp/teqp/fluiddata/dev/fluids/WATER.json");
        
        auto jig = teqp::convert_CoolProp_idealgas(FLD.dump(), 0);
        nlohmann::json ideal_gas = {{"kind","IdealHelmholtz"}, {"validate", false}, {"model", {jig}}};
        
        nlohmann::json pure = {{"components", {FLD}}};
        nlohmann::json residual = {{"kind","multifluid"}, {"validate", false}, {"model", pure}};
        auto HI = teqpflsh::properties::interfaces::teqpHelmholtzInterface(ideal_gas.dump(), residual.dump());
        
        auto getterfactory = [&HI](char k){
            return [&HI, k](double T, double rho){
                Eigen::ArrayXd z(1); z.fill(1.0);
                auto R = HI.model_residual->get_R(z);
                auto [A00, A10, A01] = HI.get_A00A10A01(T, rho, z);
                auto Ar01 = HI.model_residual->get_Ar01(T, rho, z);
                switch(k){
                    case 'H': return R*T*(1+Ar01+A10);
                    case 'S': return R*(A10-A00);
                    case 'U': return R*T*(A10);
                }
            };
        };
    
        sa.add_variable('H', getterfactory('H'));
        sa.add_variable('S', getterfactory('S'));
        sa.add_variable('U', getterfactory('U'));
        
        double Tt = 273.16;
        double Tc = 647.096;
        double eps = 1e-6;
        points = build_points(Eigen::ArrayXd::LinSpaced(300, Tt, Tc-eps), Eigen::ArrayXd::LinSpaced(300, eps, 1-eps));
    }
    void check_pair(PropertyPairs ppair){
        auto chars = pair2chars.at(ppair);
        for (auto& pt: points){
            
            // Extract the values from the point and
            // do the flash calculation
            auto val1 = pt.get_val(std::get<0>(chars));
            auto val2 = pt.get_val(std::get<1>(chars));
            auto o = sa.flash(ppair, val1, val2);
            
            CAPTURE(val1);
            CAPTURE(val2);
            CAPTURE(std::get<0>(chars));
            CAPTURE(std::get<1>(chars));
            CHECK(o);
            if (o){
                CHECK_THAT(o.value().T, WithinRel(pt.get_val('T'), 1e-4));
                CHECK_THAT(o.value().q, WithinAbs(pt.get_val('Q'), 1e-6));
            }
        }
    }
    void check_all(){
        for (auto kv: pair2chars){
            if (kv.first != PropertyPairs::PT && kv.first != PropertyPairs::HU){
                check_pair(kv.first);
            }
        }
    }
};

TEST_CASE_METHOD(SuperancillaryFixture, "Check superancillary values", "[superanc]") {
    check_all();
}
