#pragma once 

#include <chrono>

namespace teqpflsh::properties{

using namespace std::chrono;

// Sorted alphabetically (like in CoolProp)
// With 6 properties, there are 6*5/2=15 pairs
enum class PropertyPairs { ST, DT, DP, DH, DS, DU, HP, PS, PU, HS, HT, TU, SU, PT, HU };
const std::map<char, bool> uses_log = {{'T',true},{'D',true},{'P',true},{'H',false},{'S',false},{'U',false}};

const std::unordered_map<PropertyPairs, std::pair<char, char>> pair2chars = {
    {PropertyPairs::ST, {'S','T'}},
    {PropertyPairs::TU, {'T','U'}},
    {PropertyPairs::DT, {'D','T'}},
    {PropertyPairs::HT, {'H','T'}},
    {PropertyPairs::DP, {'D','P'}},
    {PropertyPairs::DH, {'D','H'}},
    {PropertyPairs::DS, {'D','S'}},
    {PropertyPairs::DU, {'D','U'}},
    {PropertyPairs::HP, {'H','P'}},
    {PropertyPairs::PS, {'P','S'}},
    {PropertyPairs::PT, {'P','T'}},
    {PropertyPairs::PU, {'P','U'}},
    {PropertyPairs::SU, {'S','U'}},
    {PropertyPairs::HU, {'H','U'}},
    {PropertyPairs::HS, {'H','S'}}
};
const std::map<std::pair<char, char>, PropertyPairs> chars2pair = {
    {{'S','T'}, PropertyPairs::ST},
    {{'D','T'}, PropertyPairs::DT},
    {{'H','T'}, PropertyPairs::HT},
    {{'T','U'}, PropertyPairs::TU},
    {{'D','P'}, PropertyPairs::DP},
    {{'D','H'}, PropertyPairs::DH},
    {{'D','S'}, PropertyPairs::DS},
    {{'D','U'}, PropertyPairs::DU},
    {{'H','P'}, PropertyPairs::HP},
    {{'P','S'}, PropertyPairs::PS},
    {{'P','T'}, PropertyPairs::PT},
    {{'P','U'}, PropertyPairs::PU},
    {{'S','U'}, PropertyPairs::SU},
    {{'H','U'}, PropertyPairs::HU},
    {{'H','S'}, PropertyPairs::HS}
};

// Get the tuple of single-character char for a pair of properties
constexpr std::pair<char, char> get_property_chars(const PropertyPairs pair){
    switch(pair){
        case PropertyPairs::ST: return {'S','T'};
        case PropertyPairs::TU: return {'T','U'};
        case PropertyPairs::DT: return {'D','T'};
        case PropertyPairs::HT: return {'H','T'};
        case PropertyPairs::DP: return {'D','P'};
        case PropertyPairs::DH: return {'D','H'};
        case PropertyPairs::DS: return {'D','S'};
        case PropertyPairs::DU: return {'D','U'};
        case PropertyPairs::HP: return {'H','P'};
        case PropertyPairs::PS: return {'P','S'};
        case PropertyPairs::PT: return {'P','T'};
        case PropertyPairs::PU: return {'P','U'};
        case PropertyPairs::SU: return {'S','U'};
        case PropertyPairs::HU: return {'H','U'};
        case PropertyPairs::HS: return {'H','S'};
    }
}
std::vector<char> get_vars(const PropertyPairs pair){
    auto chars = get_property_chars(pair);
    return {std::get<0>(chars), std::get<1>(chars)};
}
constexpr bool uses_log_(char c){
    return c == 'T' || c == 'P' || c == 'D';
}
constexpr std::tuple<bool, bool> get_relative_errors(const PropertyPairs pair){
    auto chars = get_property_chars(pair);
    return {uses_log_(std::get<0>(chars)), uses_log_(std::get<1>(chars))};
}

// Get the pair of properties given the chars for the two single properties
auto get_pair_from_chars(char ch1, char ch2){
    return chars2pair.at(std::make_pair(ch1, ch2));
}
/// Given a PropertyPair, return a tuple indicating whether the variable should be on a log basis
auto get_pair_log_scaling(const PropertyPairs pair){
    const auto chars = get_property_chars(pair);
    return std::make_tuple(uses_log.at(std::get<0>(chars)), uses_log.at(std::get<1>(chars)));
}

/**
 A generic property holder which contains arrays of numerical values for all of the thermodynamic variables
 */
template<typename Vector>
struct PropertySet{
    
    // Normally one would prefer to have these variables as stack 
    // allocations and private variables,
    // but in order to ease taking references to the allocated memory,
    // it is better to use auto ptr machinery for data allocated on
    // the heap
    const std::shared_ptr<const Vector> m_T, m_rho, m_h, m_s, m_u, m_p;
    
    // Delegate to a private constructor
    template<typename Z, typename AlphaModel>
    PropertySet(const Vector& Tvec, const Vector& rhovec, const Z& z, const AlphaModel& alpha) : PropertySet(Tvec, rhovec, build_arrays(Tvec, rhovec, z, alpha)){};
    
    PropertySet(PropertySet&& other) :
        m_T(std::move(other.m_T)),
        m_rho(std::move(other.m_rho)),
        m_h(std::move(other.m_h)),
        m_s(std::move(other.m_s)),
        m_u(std::move(other.m_u)),
        m_p(std::move(other.m_p)) {
            int rrr =0;
        };
    
    PropertySet(const PropertySet&) = default;
    
    const Vector& get_T() { return *m_T.get(); }
    const Vector& get_rho() { return *m_rho.get(); }
    const Vector& get_h() { return *m_h.get(); }
    const Vector& get_s() { return *m_s.get(); }
    const Vector& get_u() { return *m_u.get(); }
    const Vector& get_p() { return *m_p.get(); }
    
    auto get_Trho(std::size_t idx) const{ return std::make_tuple((*m_T)(idx), (*m_rho)(idx));}
    
    const std::shared_ptr<const Vector>& get_array_ptr(char k) const{
        switch(k){
            case 'T': return m_T;
            case 'D': return m_rho;
            case 'P': return m_p;
            case 'H': return m_h;
            case 'S': return m_s;
            case 'U': return m_u;
            default: throw std::invalid_argument("Bad key to get_array:"+std::to_string(k));
        }
    }
    const Vector& get_array(char k) const{
        return *get_array_ptr(k).get();
    }
    auto get_arrays(PropertyPairs pair) const{
        auto keys = get_property_chars(pair);
        return std::make_tuple(get_array(std::get<0>(keys)), get_array(std::get<1>(keys)) );
    }
    auto get_array_ptrs(PropertyPairs pair) const{
        auto keys = get_property_chars(pair);
        return std::make_tuple(get_array_ptr(std::get<0>(keys)), get_array_ptr(std::get<1>(keys)) );
    }
    
private:
    
    /// An internal data structure used to store the memory used to allocate the new pointers
    struct _newarrays{
        Vector h, s, u, p;
    };
    
    template<typename Z, typename AlphaModel>
    auto build_arrays(const Vector& Tvec, const Vector& rhovec, const Z& z, const AlphaModel& alpha){
        _newarrays n;
        n.p = Vector(Tvec.size());
        n.h = Vector(Tvec.size());
        n.s = Vector(Tvec.size());
        n.u = Vector(Tvec.size());
        
        double R = alpha.get_R(z);
        
        for (auto i = 0; i < Tvec.size(); ++i){
            double T_ = Tvec[i], Trecip = 1/T_, dTrecipdT = -Trecip*Trecip;
            double rho_ = rhovec[i];
            
            auto [A00, A10, A01] = alpha.get_A00A10A01(T_, rhovec[i], z);
            
            auto alpha = [&](){ return A00; };
            auto dalphadTrecip = [&](){ return A10/Trecip; };
            auto dalphadrho = [&](){ return A01/rho_; };
            
            auto a = [&](){ return alpha()*R*T_; };
            auto dadTrecip = [&](){ return R/(Trecip*Trecip)*(Trecip*dalphadTrecip()-alpha());};
            auto dadrho = [&]() {return R/Trecip*dalphadrho();};
            
            n.p[i] = rho_*rho_*dadrho();
            n.s[i] = Trecip*Trecip*dadTrecip();
            n.u[i] = a() + Tvec[i]*n.s[i];
            n.h[i] = n.u[i] + n.p[i]/rhovec[i];
        }
        return n;
    }
    
 // The private constructor
 PropertySet(const Vector& Tvec, const Vector& rhovec, const _newarrays &n)
    :
    m_T(std::make_shared<const Vector>(Tvec)),
    m_rho(std::make_shared<const Vector>(rhovec)),
    m_h(std::make_shared<const Vector>(n.h)),
    m_s(std::make_shared<const Vector>(n.s)),
    m_u(std::make_shared<const Vector>(n.u)),
    m_p(std::make_shared<const Vector>(n.p)) {}
};


} // namespace teqpflsh::properties
