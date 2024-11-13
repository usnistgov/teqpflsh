#pragma once

#include <optional>
#include <memory>
#include <string>
#include <variant>
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <map>

namespace teqpflsh::quadtree {

/**
 In principle it can be a bit more efficient to store center and half-width of the quad, but it is a
 bad idea for this application because if the values of the edges differ by many orders
 of magnitude, there can be a catastropic loss in precision. For instance if you have 1e-10
 as the minimum value of x and 1e9 as the max, the midpoint is 0.50000000e9, and so is the half width,
 so the minval is calculated from 0.500000e9 - 0.500000e9 = 0.000 which is not right.
 */
template <typename LeafContents>
class QuadNode{
public:
    struct Children{
        std::unique_ptr<QuadNode> NE, NW, SE, SW;
    };
public:
    std::variant<LeafContents, Children> owned;
    using OwnerOptions = decltype(owned);
    
    const double m_xmin, m_xmax, m_ymin, m_ymax;
    const std::size_t depth;
    
public:
    QuadNode(OwnerOptions&& contents, double xmin, double xmax, double ymin, double ymax, std::size_t depth)
    : owned(std::move(contents)), m_xmin(xmin), m_xmax(xmax), m_ymin(ymin), m_ymax(ymax), depth(depth) {};
    
    constexpr bool is_terminal() const { return std::holds_alternative<LeafContents>(owned); }
    
    constexpr double xmin() const { return m_xmin; }
    constexpr double xmax() const { return m_xmax; }
    constexpr double ymin() const { return m_ymin; }
    constexpr double ymax() const { return m_ymax; }
    constexpr double xmid() const { return (m_xmin + m_xmax)/2; }
    constexpr double ymid() const { return (m_ymin + m_ymax)/2; }
    
    constexpr QuadNode& NW() const { return *std::get<Children>(owned).NW; }
    constexpr QuadNode& NE() const { return *std::get<Children>(owned).NE; }
    constexpr QuadNode& SW() const { return *std::get<Children>(owned).SW; }
    constexpr QuadNode& SE() const { return *std::get<Children>(owned).SE; }
    
    /// Return true if a point is within the envelope (the bounding box)
    bool in_envelope(double x, double y) const {
        if (x < xmin() || x > xmax()){
            return false;
        }
        if (y < ymin() || y > ymax()){
            return false;
        }
        return true;
    }
    
    /**
     \param check_bounds If true, check whether an input is within the outer bounding box, will be set to false for recursive calls descending the tree because they are assumed to be valid
     */
    const QuadNode& getNode(double x, double y, bool check_bounds=true) const {
        if (!is_terminal()){
            auto& ch = std::get<Children>(owned);
            if (x > xmid()){
                if (y > ymid()){
                    return ch.NE->getNode(x, y, false);
                }
                else{
                    return ch.SE->getNode(x, y, false);
                }
            }
            else{
                if (y > ymid()){
                    return ch.NW->getNode(x, y, false);
                }
                else{
                    return ch.SW->getNode(x, y, false);
                }
            }
        }
        if (check_bounds){
            if (x < xmin() || x > xmax()){
                std::stringstream s; s << std::setprecision(17) << "x of " << x << " is out of range (" << xmin() << "," << xmax() << ")";
                throw std::invalid_argument(s.str());
            }
            if (y < ymin() || y > ymax()){
                std::stringstream s; s << std::setprecision(17) << "y of " << y << " is out of range (" << ymin() << "," << ymax() << ")";
                throw std::invalid_argument(s.str());
            }
        }
        return *this;
    }
    
    const LeafContents& get_contents() const{
        return std::get<LeafContents>(owned);
    }
    
    auto get_bisected_children(){
        Children ch;
        ch.NW = std::make_unique<QuadNode<LeafContents>>(LeafContents{}, xmin(), xmid(), ymid(), ymax(), depth+1);
        ch.SW = std::make_unique<QuadNode<LeafContents>>(LeafContents{}, xmin(), xmid(), ymin(), ymid(), depth+1);
        ch.NE = std::make_unique<QuadNode<LeafContents>>(LeafContents{}, xmid(), xmax(), ymid(), ymax(), depth+1);
        ch.SE = std::make_unique<QuadNode<LeafContents>>(LeafContents{}, xmid(), xmax(), ymin(), ymid(), depth+1);
        return ch;
    }
    
    auto bisect(){
        owned = get_bisected_children();
    }
    
    template<typename Callable>
    void walk_leavesT(const Callable& f) const{
        if (is_terminal()){
            f(*this);
        }
        else{
            auto& ch = std::get<Children>(owned);
            ch.NW->walk_leaves(f);
            ch.NE->walk_leaves(f);
            ch.SW->walk_leaves(f);
            ch.SE->walk_leaves(f);
        }
    }
    
    void walk_leaves(const std::function<void(QuadNode<LeafContents>&)>& f){
        if (is_terminal()){
            f(*this);
        }
        else{
            auto& ch = std::get<Children>(owned);
            ch.NW->walk_leaves(f);
            ch.NE->walk_leaves(f);
            ch.SW->walk_leaves(f);
            ch.SE->walk_leaves(f);
        }
    }
    void walk_leaves_const(const std::function<void(const QuadNode<LeafContents>&)>& f) const{
        if (is_terminal()){
            f(*this);
        }
        else{
            auto& ch = std::get<Children>(owned);
            ch.NW->walk_leaves_const(f);
            ch.NE->walk_leaves_const(f);
            ch.SW->walk_leaves_const(f);
            ch.SE->walk_leaves_const(f);
        }
    }
    void walk_all(const std::function<void(QuadNode<LeafContents>&)>& f){
        f(*this);
        if (!is_terminal()){
            auto& ch = std::get<Children>(owned);
            ch.NW->walk_all(f);
            ch.NE->walk_all(f);
            ch.SW->walk_all(f);
            ch.SE->walk_all(f);
        }
    }
    
    void recursive_conditional_splitter(const std::function<std::optional<Children>(QuadNode<LeafContents>&)>& f){
        if (is_terminal()){
            auto maybe_children = f(*this);
            if (maybe_children){
                owned.template emplace<Children>(std::move(maybe_children.value()));
            }
        }
        else{
            auto& ch = std::get<Children>(owned);
            ch.NW->recursive_conditional_splitter(f);
            ch.NE->recursive_conditional_splitter(f);
            ch.SW->recursive_conditional_splitter(f);
            ch.SE->recursive_conditional_splitter(f);
        }
    }
    
    void bisect_all(){
        recursive_walker([](QuadNode<LeafContents>&node){
            node.bisect();
        });
    }
    
    std::size_t count_leaves() const {
        std::size_t leaf_counter = 0;
        walk_leaves_const([&](const QuadNode<LeafContents>&node){
            if (node.is_terminal()){ leaf_counter++; };
        });
        return leaf_counter;
    }
};

}/* namespace teqpflsh::quadtree */
