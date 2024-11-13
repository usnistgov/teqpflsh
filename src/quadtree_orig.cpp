
#include <optional>
#include <memory>
#include <string>
#include <variant>
#include <iostream>
#include <functional>

class HolderThingy{
public:
    std::variant<int, std::string, std::shared_ptr<HolderThingy>> h;
    HolderThingy(const decltype(h)& h) : h(std::move(h)) {};
};

//#include <Eigen/Dense>
//
//// This uses the ray-casting algorithm to decide whether the point is inside
//// the given polygon. See https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
//bool pnpoly(const Eigen::MatrixX2d &poly, double x, double y)
//{
//    // If we never cross any lines we're inside.
//    bool inside = false;
//
//    // Loop through all the edges.
//    for (int i = 0; i < poly.rows(); ++i)
//    {
//        // i is the index of the first vertex, j is the next one.
//        // The original code uses a too-clever trick for this.
//        int j = (i + 1) % poly.rows();
//
//        // The vertices of the edge we are checking.
//        double xp0 = poly(i, 0);
//        double yp0 = poly(i, 1);
//        double xp1 = poly(j, 0);
//        double yp1 = poly(j, 1);
//
//        // Check whether the edge intersects a line from (-inf,y) to (x,y).
//
//        // First check if the line crosses the horizontal line at y in either direction.
//        if (((yp0 <= y) && (yp1 > y)) || ((yp1 <= y) && (yp0 > y)))
//        {
//            // If so, get the point where it crosses that line. This is a simple solution
//            // to a linear equation. Note that we can't get a division by zero here -
//            // if yp1 == yp0 then the above if will be false.
//            double cross = (xp1 - xp0) * (y - yp0) / (yp1 - yp0) + xp0;
//
//            // Finally check if it crosses to the left of our test point. You could equally
//            // do right and it should give the same result.
//            if (cross < x)
//                inside = !inside;
//        }
//    }
//    return inside;
//}



template <typename T>
class QuadNode{
public:
    struct Children{
        std::shared_ptr<QuadNode> NE, NW, SE, SW;
    };
public:
    std::variant<T, Children> owned;
    using OwnerOptions = decltype(owned);
    
    const double xmin, xmax, ymin, ymax;
    const std::size_t depth;
    const double xmid, ymid;
    
    bool _is_terminal;
    
public:
    QuadNode(const OwnerOptions&& contents, double xmin, double xmax, double ymin, double ymax, std::size_t depth)
    : owned(std::move(contents)), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), depth(depth), xmid((xmin+xmax)/2), ymid((ymin+ymax)/2), _is_terminal(std::holds_alternative<T>(owned)) {};
    
    bool is_terminal() const { return _is_terminal; }
    void set_terminal (bool terminal){ _is_terminal = terminal; }
    
    QuadNode& getNode(double x, double y) {
        if (!_is_terminal){
            auto& ch = std::get<Children>(owned);
            if (x > xmid){
                if (y > ymid){
                    return ch.NE->getNode(x, y);
                }
                else{
                    return ch.SE->getNode(x, y);
                }
            }
            else{
                if (y > ymid){
                    return ch.NW->getNode(x, y);
                }
                else{
                    return ch.SW->getNode(x, y);
                }
            }
        }
        return *this;
    }
    
    auto get_bisected_children(){
        Children ch;
        ch.NW = std::make_shared<QuadNode<T>>(T{}, xmin, xmid, ymid, ymax, depth+1);
        ch.SW = std::make_shared<QuadNode<T>>(T{}, xmin, xmid, ymin, ymid, depth+1);
        ch.NE = std::make_shared<QuadNode<T>>(T{}, xmid, xmax, ymid, ymax, depth+1);
        ch.SE = std::make_shared<QuadNode<T>>(T{}, xmid, xmax, ymin, ymid, depth+1);
        return ch;
    }
    
    auto bisect(){
        owned = get_bisected_children(); _is_terminal=false;
    }
    
    void recursive_walker(const std::function<void(QuadNode<T>&)>& f){
        if (is_terminal()){
            f(*this);
        }
        else{
            auto& ch = std::get<Children>(owned);
            ch.NW->recursive_walker(f);
            ch.NE->recursive_walker(f);
            ch.SW->recursive_walker(f);
            ch.SE->recursive_walker(f);
        }
    }
    
    void bisect_all(){
        recursive_walker([](QuadNode<T>&node){
            node.bisect();
        });
    }
    
    std::size_t count_leaves(){
        std::size_t leaf_counter = 0;
        recursive_walker([&](QuadNode<T>&node){
            if (node.is_terminal()){ leaf_counter++; };
        });
        return leaf_counter;
    }
};

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>

TEST_CASE("quadtree", "[quadtree]"){
    
    class NodeContents{
    public:
        auto operator ()(double x, double y){
            return x+y;
        }
    };
    
    NodeContents con;
    HolderThingy ht{2};
    
    QuadNode<NodeContents> node(con, 1, 2, 3, 4, 0);
//    SpaceFillingQuadTree<NodeContents> tree(node);
//    node.bisected_children();
    REQUIRE(node.is_terminal());
    
    QuadNode<NodeContents> next(node.get_bisected_children(), 1, 2, 3, 4, 0);
    next.bisect();
    REQUIRE(!next.is_terminal());
    
    std::function<void(QuadNode<NodeContents>&)> f = [](QuadNode<NodeContents>&node){
        std::cout << "(" << node.xmin << "," << node.xmax << ") x (" << node.ymin << "," << node.ymax << ")" << node.depth << std::endl;
    };
    next.recursive_walker(f);
    
    auto get_splitted_tree = [](int n){
        QuadNode<NodeContents> node(NodeContents{}, 1, 2, 3, 4, 0);
        for (auto split = 0; split < n; ++split){
            node.bisect_all();
        }
        return node;
    };
    auto node5 = get_splitted_tree(5);
    auto node7 = get_splitted_tree(7);
    auto node9 = get_splitted_tree(9);
    auto node11 = get_splitted_tree(11);
    
    std::cout << node5.count_leaves() << std::endl;
    std::cout << node7.count_leaves() << std::endl;
    std::cout << node9.count_leaves() << std::endl;
    
    std::cout << node7.getNode(1.7, 3.7).xmin << std::endl;
    
    next.recursive_walker(f);
    
    BENCHMARK("build tree 5"){
        return get_splitted_tree(5);
    };
    BENCHMARK("build tree 7"){
        return get_splitted_tree(7);
    };
    BENCHMARK("build tree 9"){
        return get_splitted_tree(9);
    };
    
    
    BENCHMARK("time getting node"){
        return node.getNode(1.7, 3.7).xmin;
    };
    BENCHMARK("time getting node5"){
        return node5.getNode(1.7, 3.7).xmin;
    };
    BENCHMARK("time getting node7"){
        return node7.getNode(1.7, 3.7).xmin;
    };
    BENCHMARK("time getting node9"){
        return node9.getNode(1.7, 3.7).xmin;
    };
    BENCHMARK("time getting node11"){
        return node11.getNode(1.7, 3.7).xmin;
    };
};
//
//TEST_CASE("pnpoly", "[pnpoly]"){
//    int N = 10;
//    Eigen::MatrixX2d points(N, 2);
//    auto t = Eigen::ArrayXd::LinSpaced(N, 0, 2*EIGEN_PI);
//    points.col(0) = cos(t);
//    points.col(1) = sin(t);
//    BENCHMARK("call pnpoly"){
//        return pnpoly(points, 1.7, 0.5);
//    };
//};
