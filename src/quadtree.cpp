#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "teqpflsh/polyquadtree.hpp"
#include "teqpflsh/flasher.hpp"

using namespace teqpflsh;

#include<mach/mach.h>

#include "boost/math/tools/toms748_solve.hpp"

TEST_CASE("Rootfinding", "[TOMS748]"){
    using namespace boost::math::tools;
    auto f = [](double x){ return sqrt(x)-0.7; };
    std::size_t max_iter = 100;
    auto [l,r] = toms748_solve(f, 0.0, 1.0, eps_tolerance<double>(1), max_iter);
    CHECK_THAT(r, Catch::Matchers::WithinRelMatcher(0.49, 1e-13));
}

#include "teqp/json_tools.hpp"
#include "teqpflsh/superancillary/superancillary.hpp"

TEST_CASE("Superancillary construction", "[superanc]"){
    nlohmann::json j = teqp::load_a_JSON_file("../src/testdata/WATER_exps.json");
    teqpflsh::superancillary::SuperAncillary<Eigen::ArrayXd> sa{j};
}

auto get_resident_memory = [](){
    
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

const double MYPI = 3.141592653589793238462643383279502884L;

auto create_circle_xy(int N, double x0, double y0, double r){
    std::vector<double> x, y;
    for (auto t = 0.0; t < 2*MYPI; t += 2*MYPI/(N-1)){
        x.push_back(x0+r*cos(t)); y.push_back(y0+r*sin(t));
    }
    return std::make_tuple(x, y);
};

auto create_circle(GeometryFactory::Ptr& factory, int N, double x0, double y0, double r){
    CoordinateSequence seq;
    
    for (auto t = 0.0; t < 2*MYPI; t += 2*MYPI/(N-1)){
        seq.add(x0+r*cos(t), y0+r*sin(t));
    }
    seq.closeRing();
    return std::unique_ptr<Geometry>(factory->createPolygon(std::move(seq)));
};

TEST_CASE("region2d"){
    auto [x, y] = create_circle_xy(1000, 0, 0, 1);
    QuadRegion2D reg(x, y);
    reg.do_splits(9);
    double Areainside = 0.0;
    double Areaintersection = 0.0;
    int Ninside = 0.0;
    int Nintersection = 0.0;
    auto& qtree = reg.get_quadtree_ro();
    auto& qtree_rw = reg.get_quadtree_rw();
    for (auto& leaf: qtree.get_leaves()){
        auto& c = leaf->get_contents();
        if (c.status == PQTStatus::inside){
            Areainside += c.poly->getArea();
            Ninside += 1;
        }
        if (c.status == PQTStatus::intersection){
            Areaintersection += c.poly->getArea();
            Nintersection += 1;
        }
    }
    double area_total = Areainside+Areaintersection;
    CHECK_THAT(area_total, Catch::Matchers::WithinRel(MYPI, 1e-1));
    
    BENCHMARK("contains"){
        return qtree_rw.getNode(0.9999999999, 0).get_contents().status == PQTStatus::inside;
    };
}

#include <boost/multiprecision/cpp_bin_float.hpp>
using namespace boost::multiprecision;
using my_float_type = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30U>>;

TEST_CASE("Matrix solving in higher precision"){
    
    Eigen::Matrix2d A; A << 70.0518,34.9632, 0.270827,-2.40215e-05;
    Eigen::Vector2d b; b << -5.61198e-05, 1.65601e-12;
    BENCHMARK("double precision"){
        return A.fullPivLu().solve((-b).matrix());
    };
    
    auto Amp = A.cast<my_float_type>();
    auto bmp = b.cast<my_float_type>();
    
    BENCHMARK("extended precision"){
        return Amp.fullPivLu().solve((-bmp).matrix());
    };
}


#include "geos/geom/prep/PreparedPolygon.h"
TEST_CASE("Test prepared"){
    GeometryFactory::Ptr factory = GeometryFactory::create();
    auto c = create_circle(factory, 100000, 0, 0, 1);
    geos::geom::prep::PreparedPolygon prep(c.get());
    BENCHMARK("preparation"){
        return geos::geom::prep::PreparedPolygon(c.get());
    };
    auto pt = factory->createPoint(CoordinateXY{0.00001, 0.00001});
    
    BENCHMARK("normal contains"){
        return c->contains(pt.get());
    };
    BENCHMARK("prepared contains"){
        return prep.contains(pt.get());
    };
}
TEST_CASE("pairs and keys"){
    SECTION("HS"){
        auto pair = get_pair_from_chars('H','S');
        char c1, c2;
        std::tie(c1, c2) = get_property_chars(pair);
        CHECK(c1 == 'H');
        CHECK(c2 == 'S');
    }
    SECTION("DH"){
        auto pair = get_pair_from_chars('D','H');
        char c1, c2;
        std::tie(c1, c2) = get_property_chars(pair);
        CHECK(c1 == 'D');
        CHECK(c2 == 'H');
    }
}



TEST_CASE("quadtree", "[quadtree]"){
    
    class NodeContents{
    public:
        std::unique_ptr<Geometry> poly = nullptr;
    };
    
    NodeContents con;
    
    QuadNode<NodeContents> node(std::move(con), 1, 2, 3, 4, 0);
    REQUIRE(node.is_terminal());
    
    QuadNode<NodeContents> next(node.get_bisected_children(), 1, 2, 3, 4, 0);
    next.bisect();
    REQUIRE(!next.is_terminal());
    
    std::function<void(QuadNode<NodeContents>&)> f = [](QuadNode<NodeContents>&node){
        std::cout << "(" << node.xmin() << "," << node.xmax() << ") x (" << node.ymin() << "," << node.ymax() << ")" << node.depth << std::endl;
    };
    next.walk_leaves(f);
    
    GeometryFactory::Ptr factory = GeometryFactory::create();
    
//    next.recursive_conditional_splitter(poly_splitter);
    
    auto area_stats = [](QuadNode<NodeContents>&node) -> void{
        std::vector<double> areas;
        std::function<void(QuadNode<NodeContents>&)> areasummer = [&areas](QuadNode<NodeContents>&node){
            if (node.get_contents().poly){
                areas.push_back(node.get_contents().poly->getArea());
            }
        };
        node.walk_leaves(areasummer);
        auto total = std::accumulate(areas.begin(), areas.end(), 0.0);
        
        std::cout << areas.size() << " polygons summing to " << total << ". #leaves: " << node.count_leaves() << std::endl;
    };
    
    auto base_bytes = [](QuadNode<NodeContents>&node){
        std::vector<double> bytes;
        std::function<void(QuadNode<NodeContents>&)> summer = [&bytes](QuadNode<NodeContents>&node){
            bytes.push_back(sizeof(node));
        };
        node.walk_all(summer);
        auto total = std::accumulate(bytes.begin(), bytes.end(), 0.0);
        std::cout << bytes.size() << " nodes(all) summing to " << total/1024/1024 << " MiB "<< std::endl;
    };
    
    auto poly_splitter = [&factory](QuadNode<NodeContents>&node) -> std::optional<std::decay_t<decltype(node)>::Children> {
        using node_t = std::decay_t<decltype(node)>;
        node_t::Children ch;
        
        // Do polygon intersection between the new rectangles and the polygon associated with the base node
        auto maker = [&node](double xmin, double xmax, double ymin, double ymax, std::size_t depth){
            std::unique_ptr<Geometry> newpoly;
            const NodeContents& contents = node.get_contents();
            if (contents.poly){
                auto rect = geos::operation::intersection::Rectangle(xmin, ymin, xmax, ymax);
                auto clipped = geos::operation::intersection::RectangleIntersection::clip(*(contents.poly), rect);
                auto N = clipped->getNumGeometries();
                if (clipped->getArea() > 0){
                    newpoly = clipped->getGeometryN(0)->clone();
                }
            }
            return std::make_unique<node_t>(NodeContents{std::move(newpoly)}, xmin, xmax, ymin, ymax, depth);
        };
        ch.NW = maker(node.xmin(), node.xmid(), node.ymid(), node.ymax(), node.depth+1);
        ch.SW = maker(node.xmin(), node.xmid(), node.ymin(), node.ymid(), node.depth+1);
        ch.NE = maker(node.xmid(), node.xmax(), node.ymid(), node.ymax(), node.depth+1);
        ch.SE = maker(node.xmid(), node.xmax(), node.ymin(), node.ymid(), node.depth+1);
        
        return ch;
    };
    
    auto get_splitted_tree = [&poly_splitter, &factory](int n){
        
        QuadNode<NodeContents> node(NodeContents{create_circle(factory, 1000, 0, 0, 1)}, -1, 1, -1, 1, 0);
//        QuadNode<NodeContents> node(NodeContents{}, -1, 1, -1, 1, 0);

        for (auto split = 0; split < n; ++split){
            node.recursive_conditional_splitter(poly_splitter);
        }
        return node;
    };
    
    auto r0 = get_resident_memory();
    auto [x, y] = create_circle_xy(1000, 0, 0, 1);
    PolyQuadTree polytree(x,y);
    polytree.do_splits(3);
    polytree.area_stats();
    auto r1 = get_resident_memory();
    std::cout << std::setprecision(14) << double(r1-r0)/1024/1024 << " MiB per approximated circle" << std::endl;
    
    BENCHMARK("makeit"){
        auto [x, y] = create_circle_xy(1000, 0, 0, 1);
        PolyQuadTree polytree(x,y);
        polytree.do_splits(2);
        return polytree;
    };
    
    auto resident0 = get_resident_memory();
    int Npts = 300;
    std::vector<std::unique_ptr<Geometry>> rr;
    for (auto i = 0U; i < 10000; ++i){
        rr.emplace_back(create_circle(factory, Npts, 0, 0, 1));
    }
    auto resident1 = get_resident_memory();
    std::cout << double(resident1-resident0)/10000 << " bytes per circle" << std::endl;
    std::cout << double(resident1-resident0)/10000-(Npts+1)*16 << " bytes overhead per circle" << std::endl;
    
    using QT = QuadNode<NodeContents>;
    std::cout << sizeof(NodeContents) << std::endl;
    std::cout << sizeof(QT::Children) << std::endl;
    std::cout << sizeof(std::variant<NodeContents, QT::Children>) << std::endl;
    std::cout << sizeof(QT) << std::endl;
    
    std::unordered_map<std::size_t, QuadNode<NodeContents>> trees;
    for (std::size_t nsplit = 0U; nsplit < 10; nsplit += 1){
        auto resident0 = get_resident_memory();
        QuadNode<NodeContents> tree = get_splitted_tree(nsplit);
        area_stats(tree);
        base_bytes(tree);
        auto resident1 = get_resident_memory();
        std::cout << std::setprecision(12) << double(resident1-resident0)/1024/1024 << " MiB" << std::endl;
        
//        trees.insert({nsplit, std::move(tree)});
    }
    
    BENCHMARK("build tree 5"){
        return get_splitted_tree(5);
    };
    BENCHMARK("build tree 7"){
        return get_splitted_tree(7);
    };
    BENCHMARK("build tree 9"){
        return get_splitted_tree(9);
    };
    
    
    
//    BENCHMARK("time getting node"){
//        return trees[0].getNode(1.7, 3.7).xmin();
//    };
//    BENCHMARK("time getting node5"){
//        return trees[5].getNode(1.7, 3.7).xmin();
//    };
//    BENCHMARK("time getting node7"){
//        return trees[7].getNode(1.7, 3.7).xmin();
//    };
//    BENCHMARK("time getting node9"){
//        return trees[9].getNode(1.7, 3.7).xmin();
//    };
//    BENCHMARK("time getting node11"){
//        return trees[11].getNode(1.7, 3.7).xmin();
//    };
};
