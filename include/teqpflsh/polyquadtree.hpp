#pragma once

#include "teqpflsh/trees/quadtree/quadtree.hpp"

/* For geometry operations */
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/intersection/Rectangle.h>
#include <geos/operation/intersection/RectangleIntersection.h>

namespace teqpflsh {

/* Geometry/GeometryFactory */
using namespace geos::geom;

enum class PQTStatus {inside, outside, intersection};

class PolyQuadTree{
    
public:
    class LeafContents{
    public:
        std::unique_ptr<Geometry> poly = nullptr;
        PQTStatus status = PQTStatus::intersection;
    };
    
    using QuadTree = quadtree::QuadNode<LeafContents>;
    
private:
    auto create_base_polygon(const std::vector<double>&x, const std::vector<double>&y){
        CoordinateSequence seq;
        for (auto i=0; i < x.size(); ++i){
            seq.add(x[i], y[i]);
        }
        seq.closeRing();
        return std::unique_ptr<Geometry>(factory->createPolygon(std::move(seq)));
    }
    auto create_tree(const std::vector<double>&x, const std::vector<double>&y){
        auto base_poly = create_base_polygon(x, y);
        
        auto coords = base_poly->getEnvelope()->getCoordinates();
//        static_assert(coords->size() == 5U);
        
        // 5 points, starting in lower left and going counter-clockwise
        // so max is at index 2
        auto xmin = coords->getX(0);
        auto ymin = coords->getY(0);
        auto xmax = coords->getX(2);
        auto ymax = coords->getY(2);
        
        return QuadTree(LeafContents{std::move(base_poly)}, xmin, xmax, ymin, ymax, 0);
    }
    auto create_tree_poly(std::unique_ptr<CoordinateSequence>&& seq){
        
        std::unique_ptr<geos::geom::LinearRing> lr = factory->createLinearRing(std::move(seq));
        std::unique_ptr<Geometry> newpoly = factory->createPolygon(std::move(lr));
        
        auto coords = newpoly->getEnvelope()->getCoordinates();
//        static_assert(coords->size() == 5U);
        
        // 5 points, starting in lower left and going counter-clockwise
        // so max is at index 2
        auto xmin = coords->getX(0);
        auto ymin = coords->getY(0);
        auto xmax = coords->getX(2);
        auto ymax = coords->getY(2);
        
        return QuadTree{LeafContents{std::move(newpoly)}, xmin, xmax, ymin, ymax, 0};
    }
public:
    
    GeometryFactory::Ptr factory;
    QuadTree tree;
    
    PolyQuadTree(const std::vector<double>&x, const std::vector<double>&y) : factory(GeometryFactory::create()), tree(create_tree(x, y)){}
    PolyQuadTree(std::unique_ptr<CoordinateSequence>&& seq): factory(GeometryFactory::create()), tree(create_tree_poly(std::move(seq))){}
    
    /// Return true if a point is within the envelope (the bounding box)
    bool in_envelope(double x, double y) const {
        return tree.in_envelope(x, y);
    }
    /// Return true if a point is inside a node that is fully inside the polygon or intersects the polygon
    bool is_inside_or_inter(double x, double y) const {
        const auto& node = tree.getNode(x, y, true /* check_bounds*/);
        auto status = node.get_contents().status;
        return status == PQTStatus::inside || status == PQTStatus::intersection;
    }
    /// Get a node of the tree
    const auto& getNode(double x, double y, bool check_bounds=true) const {
        return tree.getNode(x, y, check_bounds);
    }
    /// Print to stdout some information about the nodes of the tree
    auto area_stats(){
        std::vector<double> areas;
        std::function<void(QuadTree&)> areasummer = [&areas](QuadTree&node){
            if (node.get_contents().poly){
                areas.push_back(node.get_contents().poly->getArea());
            }
        };
        tree.walk_leaves(areasummer);
        auto total = std::accumulate(areas.begin(), areas.end(), 0.0);
        std::cout << areas.size() << " polygons summing to " << total << ". #leaves: " << tree.count_leaves() << std::endl;
    };
    
    auto do_splits(int Nsplits){
        auto poly_splitter = [this](quadtree::QuadNode<LeafContents>&node) -> std::optional<std::decay_t<decltype(node)>::Children> {
            const auto& poly = node.get_contents().poly;
            if (!poly || node.get_contents().status == PQTStatus::inside || node.get_contents().status == PQTStatus::outside){ return std::nullopt; }
            
            using node_t = std::decay_t<decltype(node)>;
            
            auto rect = geos::operation::intersection::Rectangle(node.xmin(), node.ymin(), node.xmax(), node.ymax());
            auto clipped = geos::operation::intersection::RectangleIntersection::clip(*poly, rect);
            
            // If intersection area is greater than zero, there is some intersection
            // and either the rectangle is entirely within the polygon
            // or it intersects the polygon
            //
            //
            
            if (clipped->getArea() > 0.0){
                
                // Do polygon intersection between the new rectangles and the polygon associated with the base node
                auto make_child = [&poly, this](double xmin, double xmax, double ymin, double ymax, std::size_t depth){
                    auto rect = geos::operation::intersection::Rectangle(xmin, ymin, xmax, ymax);
                    auto clipped = geos::operation::intersection::RectangleIntersection::clip(*poly, rect);
//                    auto clipped = poly->intersection(rect.toPolygon(*factory.get()).get());
                    auto Arect = rect.toPolygon(*factory.get())->getArea();
                    auto Aclipped = clipped->getArea();
                    
                    auto iswithin = poly->contains(rect.toPolygon(*factory.get()).get());
                    
                    // There can be degeneracy when you have a rectangle that is aligned with an edge of the
                    // polygon, such that although it is actually entirely within, rounding takes some points
                    // outside the envelope, so make the within check based upon a tolerance on the area.
                    // if the area is identical, then you are done
                    
//                    bool iscomplete = iswithin || std::abs(Aclipped-Arect) < 1e-8*Arect;
                    PQTStatus status = PQTStatus::intersection;
                    if (iswithin){
                        status = PQTStatus::inside;
                    }
                    else if (Aclipped == 0){
                        status = PQTStatus::outside;
                    }
                    std::unique_ptr<Geometry> newpoly; // No new poly
                    if (Aclipped != 0.0){
                        newpoly.swap(clipped);
                    }
                    return std::make_unique<node_t>(LeafContents{std::move(newpoly), status}, xmin, xmax, ymin, ymax, depth);
                };
                node_t::Children ch;
                ch.NW = make_child(node.xmin(), node.xmid(), node.ymid(), node.ymax(), node.depth+1);
                ch.SW = make_child(node.xmin(), node.xmid(), node.ymin(), node.ymid(), node.depth+1);
                ch.NE = make_child(node.xmid(), node.xmax(), node.ymid(), node.ymax(), node.depth+1);
                ch.SE = make_child(node.xmid(), node.xmax(), node.ymin(), node.ymid(), node.depth+1);
                return ch;
            }
            return std::nullopt;
        };
        for (auto i = 0; i < Nsplits; ++i){
            tree.recursive_conditional_splitter(poly_splitter);
        }
    }
    ///
    auto get_polygon_xy(QuadTree&node) -> std::optional<std::tuple<std::vector<double>,std::vector<double>>>{
        const LeafContents& contents = node.get_contents();
        if (contents.poly){
            std::vector<CoordinateXY> coordsxy;
            contents.poly->getCoordinates()->toVector(coordsxy);
            std::vector<double> x, y;
            for (auto&&xy : coordsxy){
                x.push_back(xy.x);
                y.push_back(xy.y);
            }
            return std::make_tuple(x, y);
        }
        return std::nullopt;
    }
    auto is_complete(const QuadTree&node) const{
        return node.get_contents().status == PQTStatus::inside;
    }
    auto is_intersection(const QuadTree&node) const{
        return node.get_contents().status == PQTStatus::intersection;
    }
    auto get_status(const QuadTree&node) const{
        return node.get_contents().status;
    }
    auto get_leaves() const{
        std::vector<const QuadTree*> leaves;
        std::function<void(const QuadTree&)> walker = [&leaves](const QuadTree&node){
            leaves.push_back(&node);
        };
        tree.walk_leaves_const(walker);
        return leaves;
    }
    auto walk_leaves_const(const std::function<void(const QuadTree&)>& f) const{
        std::cout << "Into main call" << std::endl;
        tree.walk_leaves_const(f);
    }
    template<typename Callable>
    auto walk_leavesT(const Callable& f){
        std::cout << "Into main call" << std::endl;
        tree.walk_leavesT<Callable>(f);
    }
};

}/* */
