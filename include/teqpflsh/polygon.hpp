#pragma once

#include "teqpflsh/polymath.hpp"

namespace teqpflsh{

template<typename ArrayType>
class Polygon2D{
public:
	ArrayType m_x, m_y; // ArrayType must have .size() and .data() methods (so Eigen types or nanobind objects both work, but not STL types)
	auto x() const { return m_x; }
	auto y() const { return m_y; }
	bool inpoly(double x, double y){
		return pnpoly(m_x.size(), m_x.data(), m_y.data(), x, y);
	}
    
	// template<RetType>
	// RetType inpolyvec(const ArrayType& x, const ArrayType& y){
	// 	RetType out; out.resize()
	// 	return pnpoly(m_x.size(), m_x.data(), m_y.data(), x, y);
	// }
};

}
