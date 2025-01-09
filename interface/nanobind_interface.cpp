#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/ndarray.h>
#include <nanobind/eigen/dense.h>

#include "geos/geom/prep/PreparedGeometryFactory.h"
#include "geos/simplify/DouglasPeuckerSimplifier.h"
#include "geos/simplify/TopologyPreservingSimplifier.h"
#include "geos/operation/valid/MakeValid.h"

#define BOOST_MATH

#include "teqpflsh/polyquadtree.hpp"
#include "teqpflsh/flasher.hpp"
#include "teqpflsh/superancillary/superancillary.hpp"

namespace nb = nanobind;
using namespace nb::literals;

using namespace teqpflsh;

using ArrayType = teqpflsh::ArrayType;

using tensor1d = nb::ndarray<double, nb::shape<-1>, nb::c_contig, nb::device::cpu>;
using tensor1i = nb::ndarray<int, nb::shape<-1>, nb::c_contig, nb::device::cpu>;
using tensor33d = nb::ndarray<double, nb::shape<3,3>, nb::c_contig, nb::device::cpu>;

int add(int a, int b) { return a + b; }
double _indexer(const tensor1d x, std::size_t i){
    return x(i);
}

double _indexer33(const tensor33d x, std::size_t i, std::size_t j, char k){
    return x(i, j);
}

NB_MODULE(_teqpflsh_impl, m) {
    m.attr("__version__") = TEQPFLSHVERSION;
    
    using namespace properties;
    using namespace quadtree;
    using namespace kdtree;

    m.def("add", &add);
    m.def("indexer", &_indexer);
    m.def("indexer33", &_indexer33);
    
    
    //  *****************************************************************************
    //  *****************************************************************************
    //  **************************        GEOS       ********************************
    //  *****************************************************************************
    //  *****************************************************************************
    
    // Wrapping of some elements from GEOS that are handy for debugging
    
    // The GeometryFactory::Ptr of GEOS has a custom deleter that nanobind cannot wrap natively
    // so a holder class is used to control its lifetime
    struct GeometryFactoryHolder{
        GeometryFactory::Ptr factory;
        GeometryFactoryHolder() : factory(GeometryFactory::create()) {}
        auto createPolygon(CoordinateSequence&& seq){ return std::unique_ptr<Geometry>(factory->createPolygon(std::move(seq))); }
        auto createPoint(double x, double y){ return std::unique_ptr<Point>(factory->createPoint(CoordinateXY{x,y})); }
    };
    
    nb::class_<GeometryFactoryHolder>(m, "GeometryFactoryHolder")
        .def(nb::init<>())
        .def("createPolygon", &GeometryFactoryHolder::createPolygon)
        .def("createPoint", &GeometryFactoryHolder::createPoint)
        .def("makeclosedpolygon", [](GeometryFactoryHolder *holder, const tensor1d& x, const tensor1d& y){
            auto seq = CoordinateSequence();
            for (auto i = 0; i < x.size(); ++i){
                seq.add(x(i), y(i));
            }
            seq.closeRing(true);
            return holder->createPolygon(std::move(seq));
        }, "x"_a, "y"_a, "A convenience function to make a closed polygon given numpy arrays")
    ;
    
    nb::class_<GeometryFactory>(m, "GeometryFactory")
        .def("createPolygon", &GeometryFactory::createGeometry)
    ;
    
    nb::class_<CoordinateSequence>(m, "CoordinateSequence")
        .def(nb::init<>())
        .def("add", [](CoordinateSequence &cs, double x, double y) { cs.add(x,y); })
        .def("getX", &CoordinateSequence::getX)
        .def("getY", &CoordinateSequence::getY)
        .def("closeRing", &CoordinateSequence::closeRing)
        .def("getSize", &CoordinateSequence::getSize)
    ;
    
    using namespace geos::geom::prep;
    nb::class_<PreparedGeometry>(m, "PreparedGeometry")
        .def("contains", &PreparedGeometry::contains)
        .def("nearestPoints", &PreparedGeometry::nearestPoints)
    ;
    
    nb::class_<Geometry>(m, "Geometry")
        .def("getNumGeometries", &Geometry::getNumGeometries)
        .def("getGeometryN", &Geometry::getGeometryN, nb::rv_policy::reference)
        .def("intersection", &Geometry::intersection)
        .def("difference", &Geometry::difference)
        .def("DelaunayTriangulate", [](const Geometry* geo){ return geos::triangulate::polygon::ConstrainedDelaunayTriangulator::triangulate(geo); })
        .def("fastTriangulate", [](const Geometry* geo){ return geos::triangulate::polygon::PolygonTriangulator::triangulate(geo); })
    
        .def("getNumPoints", &Geometry::getNumPoints)
        .def("getXY", [](const Geometry* self){
            auto coords = self->getCoordinates();
            Eigen::ArrayXd X(coords->size()), Y(coords->size());
            for (auto j = 0; j < coords->size(); ++j){
                X(j) = coords->getX(j);
                Y(j) = coords->getY(j);
            }
            return std::make_tuple(X, Y);
        }, "Convenience function to return the X, Y coordinates as numpy arrays in Python")
        .def("containsPoint", [](const Geometry* self, const std::unique_ptr<Point>& pt){ return self->contains(pt.get()); })
        .def("getCoordinates", &Geometry::getCoordinates)
        .def_prop_ro("isValid", &Geometry::isValid)
        .def_prop_ro("isSimple", &Geometry::isSimple)
        .def("get_PreparedGeometry", [](const Geometry* self){
            return geos::geom::prep::PreparedGeometryFactory::prepare(self);
        })
        .def("getCentroid",  [](Geometry* t) { return t->getCentroid(); })
        .def("run_DouglasPeuckerSimplifier", [](const Geometry* self, double tolerance){
            geos::simplify::DouglasPeuckerSimplifier sim(self);
            sim.setDistanceTolerance(tolerance);
            return sim.getResultGeometry();
        }, "tolerance"_a)
        .def("run_TopologyPreservingSimplifier", [](const Geometry* self, double tolerance){
            geos::simplify::TopologyPreservingSimplifier sim(self);
            sim.setDistanceTolerance(tolerance);
            return sim.getResultGeometry();
        }, "tolerance"_a)
        .def("make_valid", [](const Geometry* self){
            return geos::operation::valid::MakeValid().build(self);
        })
    ;
    
    nb::class_<Point, Geometry>(m, "Point")
        .def("getX", &Point::getX)
        .def("getY", &Point::getY)
    ;
    
    //  *****************************************************************************
    //  *****************************************************************************
    //  ************************        K-D Tree       ******************************
    //  *****************************************************************************
    //  *****************************************************************************
    
    nb::class_<AbstractScaler<double>>(m, "AbstractScaler");
    nb::class_<MinMaxScaler<double>, AbstractScaler<double>>(m, "MinMaxScaler")
        .def(nb::init<double, double>(), "x_min"_a, "x_max"_a);
    nb::class_<MinMaxLogScaler<double>, AbstractScaler<double>>(m, "MinMaxLogScaler")
        .def(nb::init<double, double>(), "x_min"_a, "x_max"_a);
    
    using L2Tree_t = L2Tree<Eigen::ArrayXd>;
    nb::class_<L2Tree_t>(m, "L2Tree")
        .def("get_nearest_indexd2", nb::overload_cast<const double, const double>(&L2Tree_t::get_nearest_indexd2, nb::const_), "x"_a, "y"_a)
        .def("get_nearest_indexd2", nb::overload_cast<const EArray2&>(&L2Tree_t::get_nearest_indexd2, nb::const_), "pt"_a)
        .def("get_nearest_indexd2_many", &L2Tree_t::get_nearest_indexd2_many<tensor1d, tensor1i>, "x"_a, "y"_a, "idx"_a, "d2"_a)
        .def("get_used_bytes", &L2Tree_t::get_used_bytes)
    ;
    
    struct L2Holder{
        using num_t = double;
        const std::shared_ptr<const Eigen::ArrayXd> xptr, yptr;
        std::shared_ptr<L2Tree_t> tree;
        const std::shared_ptr<AbstractScaler<num_t>> xscaler, yscaler;
        L2Holder(const Eigen::ArrayXd&x, const Eigen::ArrayXd&y, int depth, const std::shared_ptr<AbstractScaler<double>>&xscaler, const std::shared_ptr<AbstractScaler<double>>&yscaler) 
        : 
        xptr(std::make_shared<const Eigen::ArrayXd>(x)),
        yptr(std::make_shared<const Eigen::ArrayXd>(y)),
        tree(std::make_shared<L2Tree_t>(xptr, yptr, depth, xscaler, yscaler)) {}
    };
    nb::class_<L2Holder>(m, "L2TreeHolder")
        .def(nb::init<const Eigen::ArrayXd&, const Eigen::ArrayXd&, int, const std::shared_ptr<AbstractScaler<double>>&, const std::shared_ptr<AbstractScaler<double>>&>(), "x"_a, "y"_a, "tree_depth"_a, nb::arg("xscaler") = nb::none(), nb::arg("yscaler") = nb::none() )
        .def_rw("tree", &L2Holder::tree, nb::rv_policy::reference)
        ;
    
    
    
    
    
    using PolyQuadNode = PolyQuadTree::QuadTree;
    nb::class_<PolyQuadNode>(m, "PolyQuadNode")
     .def("xmin", &PolyQuadNode::xmin)
     .def("xmax", &PolyQuadNode::xmax)
     .def("ymin", &PolyQuadNode::ymin)
     .def("ymax", &PolyQuadNode::ymax)
     .def("getNode", &PolyQuadNode::getNode, nb::rv_policy::reference)
     .def("get_contents", &PolyQuadNode::get_contents, nb::rv_policy::reference)
     .def_prop_ro("terminal", [](PolyQuadNode &t) { return t.is_terminal(); })
     .def_prop_ro("NW", [](PolyQuadNode &t) -> PolyQuadNode & { return t.NW(); } , nb::rv_policy::reference)
     .def_prop_ro("NE", [](PolyQuadNode &t) -> PolyQuadNode & { return t.NE(); } , nb::rv_policy::reference)
     .def_prop_ro("SW", [](PolyQuadNode &t) -> PolyQuadNode & { return t.SW(); } , nb::rv_policy::reference)
     .def_prop_ro("SE", [](PolyQuadNode &t) -> PolyQuadNode & { return t.SE(); } , nb::rv_policy::reference)
    ;
    
    using LeafContents = PolyQuadTree::LeafContents;
    nb::class_<LeafContents>(m, "LeafContents")
//        .def_ro("poly", &LeafContents::poly, nb::rv_policy::reference)
        .def_ro("status", &LeafContents::status)
        ;
    
    nb::class_<PolyQuadTree>(m, "PolyQuadTree")
     .def(nb::init<const std::vector<double>&, const std::vector<double>&>())
     .def("do_splits", &PolyQuadTree::do_splits)
     .def("area_stats", &PolyQuadTree::area_stats)
     .def("get_polygon_xy", &PolyQuadTree::get_polygon_xy)
     .def("is_complete", &PolyQuadTree::is_complete)
     .def("is_intersection", &PolyQuadTree::is_intersection)
    .def("get_status", &PolyQuadTree::get_status)
    .def("get_leaves", &PolyQuadTree::get_leaves, nb::rv_policy::reference)
//     .def("walk_leaves_const", &PolyQuadTree::walk_leaves_const, nb::rv_policy::reference)
//     .def("walk_leaves", &PolyQuadTree::walk_leavesT<nb::callable>)
     .def_ro("tree", &PolyQuadTree::tree)
    ;
    
    nb::enum_<PQTStatus>(m, "PQTStatus")
        .value("inside", PQTStatus::inside)
        .value("outside", PQTStatus::outside)
        .value("intersection", PQTStatus::intersection)
    ;
    
    nb::class_<QuadRegion2D>(m, "QuadRegion2D")
     .def(nb::init<const std::vector<double>, const std::vector<double>>(), nb::kw_only(), "x"_a, "y"_a)
     .def("do_splits", &QuadRegion2D::do_splits)
     .def("get_envelope", &QuadRegion2D::get_envelope)
     .def("get_quadtree_ro", &QuadRegion2D::get_quadtree_ro, nb::rv_policy::reference)
     .def("get_quadtree_rw", &QuadRegion2D::get_quadtree_rw, nb::rv_policy::reference)
     .def("do_Delaunay_triangulation", &QuadRegion2D::do_Delaunay_triangulation, nb::rv_policy::reference)
     .def("do_fast_triangulation", &QuadRegion2D::do_fast_triangulation, nb::rv_policy::reference)
     .def_prop_ro("bounding_polygon", &QuadRegion2D::get_bounding_polygon, nb::rv_policy::reference_internal)
     .def("sample_random", &QuadRegion2D::sample_random<tensor1d>)
     .def("sample_gridded", &QuadRegion2D::sample_gridded<tensor1d>)
     .def("get_coords_xy", &QuadRegion2D::get_coords_xy)
     .def("get_nonsimple_xy", &QuadRegion2D::get_nonsimple_xy)
     .def("sample_gridded_w_tree", &QuadRegion2D::sample_gridded_w_tree<tensor1d>)
    ;
    
    m.def("sample_random", &::sample_random<tensor1d>);
    
    using env = QuadRegion2D::Envelope;
    nb::class_<env>(m, "Envelope")
        .def_ro("x_min", &env::x_min)
        .def_ro("y_min", &env::y_min)
        .def_ro("x_max", &env::x_max)
        .def_ro("y_max", &env::y_max)
    ;
    
    using tr = ThermodynamicRegion<Eigen::ArrayXd>;
    nb::class_<tr>(m, "ThermodynamicRegion")
        .def_prop_ro("transformed_regions", &tr::get_transformed_regions, nb::rv_policy::reference)
        .def_prop_ro("propset_bounding", &tr::get_propset_bounding, nb::rv_policy::reference)
        .def_prop_ro("propset_Trhogrid", &tr::get_propset_Trhogrid, nb::rv_policy::reference)
        .def("add_pair", &tr::add_pair, nb::kw_only(), "proppair"_a, "Nsplit"_a, "and_kdtree"_a = true)
        .def("has_pair", &tr::has_pair)
        .def("get_transformed_region", &tr::get_transformed_region)
        .def("get_kdtree", &tr::get_kdtree)
        .def("get_starting_Trho_many", &tr::get_starting_Trho_many<tensor1d>, "proppair"_a, "val1"_a, "val2"_a, "T"_a, "rho"_a, "d2"_a)
    ;
    
    using ps = PropertySet<Eigen::ArrayXd>;
    nb::class_<ps>(m, "PropertySet")
        .def_prop_ro("T", &ps::get_T, nb::rv_policy::reference)
        .def_prop_ro("p", &ps::get_p, nb::rv_policy::reference)
        .def_prop_ro("h", &ps::get_h, nb::rv_policy::reference)
        .def_prop_ro("s", &ps::get_s, nb::rv_policy::reference)
        .def_prop_ro("u", &ps::get_u, nb::rv_policy::reference)
        .def_prop_ro("rho", &ps::get_rho, nb::rv_policy::reference)
        .def("get_array", &ps::get_array, nb::rv_policy::reference)
        .def("get_arrays", &ps::get_arrays)
    ;
    
    nb::enum_<PropertyPairs>(m, "PropertyPairs")
        .value("HS", PropertyPairs::HS)
        .value("ST", PropertyPairs::ST)
        .value("HT", PropertyPairs::HT)
        .value("TU", PropertyPairs::TU)
        .value("DP", PropertyPairs::DP)
        .value("PS", PropertyPairs::PS)
        .value("DH", PropertyPairs::DH)
        .value("DS", PropertyPairs::DS)
        .value("HP", PropertyPairs::HP)
        .value("PU", PropertyPairs::PU)
        .value("PT", PropertyPairs::PT)
        .value("SU", PropertyPairs::SU)
        .value("HU", PropertyPairs::HU)
        .value("DT", PropertyPairs::DT)
        .value("DU", PropertyPairs::DU)
    ;
    
    using scr = teqp::iteration::StoppingConditionReason;
    nb::enum_<scr>(m, "StoppingConditionReason")
        .value("fatal", scr::fatal)
        .value("success", scr::success)
        .value("keep_going", scr::keep_going);
    
    nb::class_<teqp::iteration::StoppingCondition>(m, "StoppingCondition");
    
    nb::class_<teqp::iteration::MaxAbsErrorCondition, teqp::iteration::StoppingCondition>(m, "MaxAbsErrorCondition")
     .def(nb::init<double>(), "threshold"_a)
    ;
    nb::class_<teqp::iteration::NanXDXErrorCondition, teqp::iteration::StoppingCondition>(m, "NanXDXErrorCondition")
        .def(nb::init<>())
    ;
    
    using nri = teqp::iteration::NRIterator;
    nb::class_<nri>(m, "NRIterator")
        .def_rw("verbose", &nri::verbose)
        .def("get_T", &nri::get_T)
        .def("get_rho", &nri::get_rho)
        .def("get_vals", &nri::get_vals)
        .def("reset", &nri::reset)
    
        .def("calc_step", &nri::calc_step)
        .def("calc_just_step", &nri::calc_just_step)
        .def("get_nonconstant_indices", &nri::get_nonconstant_indices)
    
        .def("calc_r", &nri::calc_r)
        .def("calc_J", &nri::calc_J)
        .def("calc_vals", &nri::calc_vals)
        .def("calc_maxabsr", &nri::calc_maxabsr)
        .def("get_maxabsr", &nri::get_maxabsr)
        
        .def("take_steps", &nri::take_steps, "N"_a, "apply_stopping"_a)
        .def("path_integration", &nri::path_integration)
     
       .def("get_step_count", &nri::get_step_count)
    ;
    
    nb::class_<RegionedFlashReturn>(m, "RegionedFlashReturn")
        .def_ro("T", &RegionedFlashReturn::T, "Temperature, K")
        .def_ro("rho", &RegionedFlashReturn::rho, "Molar density, mol/m3")
        .def_ro("reason", &RegionedFlashReturn::reason, "Enumerated value for stopping reason")
        .def_ro("step_count", &RegionedFlashReturn::step_count, "How many Newton steps were takenm")
        .def_ro("maxabsr", &RegionedFlashReturn::maxabsr, "Maximum absolute residual")
        .def_ro("msg", &RegionedFlashReturn::msg, "Message associated with stoppping reason")
        .def_ro("newton_duration_us", &RegionedFlashReturn::newton_duration_us, "How long the Newton part took, in microseconds")
        .def_ro("total_duration_us", &RegionedFlashReturn::total_duration_us, "How long the total calculation took, in microseconds")
        .def_ro("candidate_duration_us", &RegionedFlashReturn::candidate_duration_us, "How long the candidate determination part took, in microseconds")
        ;
    
    nb::class_<RegionedFlasher>(m, "RegionedFlasher")
     .def(nb::init<const std::string&, const std::string&, const ArrayType&>(), nb::kw_only(), "ideal_gas"_a, "resid"_a, "mole_fractions"_a)
     .def("add_region", &RegionedFlasher::add_region, nb::kw_only(), "T"_a, "rho"_a, "NT"_a, "Nrho"_a, "Add a region to the set of regions")
     .def("remove_all_regions", &RegionedFlasher::remove_all_regions, "Remove all the regions to restore object to its initial state")
     .def("get_regions_ro", &RegionedFlasher::get_regions_ro, nb::rv_policy::reference, "Get a read-only view of the regions")
     .def("get_regions_rw", &RegionedFlasher::get_regions_rw, nb::rv_policy::reference, "Get read-write access to the regions")
     .def("get_quadtree_intersections", &RegionedFlasher::get_quadtree_intersections, nb::rv_policy::reference)
     .def("get_NRIterator", &RegionedFlasher::get_NRIterator, "Construct a Newton iterator object")
     .def("get_starting_Trho", &RegionedFlasher::get_starting_Trho, nb::rv_policy::reference, "Get the starting temperature, density pair from the K-D tree")
     
     .def("flash", &RegionedFlasher::flash, "proppair"_a, "val1"_a, "val2"_a, "Do a flash calculation")
     .def("flash_many", &RegionedFlasher::flash_many<tensor1d>, "proppair"_a, "val1"_a, "val2"_a, "T"_a, "rho"_a, "steps"_a, "maxabs"_a, "newtontime"_a, "candtime"_a, "Do many flash calculations, for testing in Python")
    ;
    
    using ce = teqpflsh::superancillary::ChebyshevExpansion<ArrayType>;
    nb::class_<ce>(m, "ChebyshevExpansion")
     .def(nb::init<double, double, const ArrayType&>(), nb::kw_only(), "xmin"_a, "xmax"_a, "coeff"_a)
     .def_prop_ro("xmin", &ce::xmin)
     .def_prop_ro("xmax", &ce::xmax)
     .def_prop_ro("coeff", &ce::coeff, nb::rv_policy::reference)
     .def("eval", &ce::eval<double>)
     .def("eval_many", &ce::eval_many<tensor1d>)
     .def("eval_Eigen", &ce::eval_many<tensor1d>)
     .def("solve_for_x", &ce::solve_for_x)
     .def("solve_for_x_count", &ce::solve_for_x_count)
     .def("solve_for_x_many", &ce::solve_for_x_many<tensor1d>)
    ;
    
    using ca = teqpflsh::superancillary::ChebyshevApproximation1D<ArrayType>;
    nb::class_<ca>(m, "ChebyshevApproximation1D")
     .def(nb::init<std::vector<ce>&&>(), nb::kw_only(), "expansions"_a)
        .def("get_x_for_y", &ca::get_x_for_y, nb::kw_only(), "y"_a, "bits"_a = 64, "max_iter"_a=100, "boundsftol"_a=1e-13)
        .def("eval", &ca::eval)
        .def_ro("xmin", &ca::xmin)
        .def_ro("xmax", &ca::xmax)
        .def("eval_many", &ca::eval_many<tensor1d>)
        .def("count_x_for_y_many", &ca::count_x_for_y_many<tensor1d>)
        .def_prop_ro("expansions", &ca::get_expansions)
        .def_prop_ro("x_at_extrema", &ca::get_x_at_extrema)
        .def_prop_ro("monotonic_intervals", &ca::get_monotonic_intervals)
        .def("get_intervals_containing_y", &ca::get_intervals_containing_y)
    ;
    
    using teqpflsh::superancillary::MonotonicExpansionMatch;
    nb::class_<MonotonicExpansionMatch>(m, "MonotonicExpansionMatch")
        .def_ro("xmin", &MonotonicExpansionMatch::xmin)
        .def_ro("xmax", &MonotonicExpansionMatch::xmax)
        .def_ro("ymin", &MonotonicExpansionMatch::ymin)
        .def_ro("ymax", &MonotonicExpansionMatch::ymax)
        .def_ro("idx", &MonotonicExpansionMatch::idx)
    ;
    
    using teqpflsh::superancillary::IntervalMatch;
    nb::class_<IntervalMatch>(m, "IntervalMatch")
        .def_ro("xmin", &IntervalMatch::xmin)
        .def_ro("xmax", &IntervalMatch::xmax)
        .def_ro("ymin", &IntervalMatch::ymin)
        .def_ro("ymax", &IntervalMatch::ymax)
        .def_ro("expansioninfo", &IntervalMatch::expansioninfo)
    ;
    
    using satps = teqpflsh::superancillary::SuperAncillaryTwoPhaseSolution;
    nb::class_<satps>(m, "SuperAncillaryTwoPhaseSolution")
    .def_ro("T", &satps::T)
    .def_ro("q", &satps::q)
    .def_ro("counter", &satps::counter)
    ;
    
    using sa = teqpflsh::superancillary::SuperAncillary<ArrayType>;
    nb::class_<sa>(m, "SuperAncillary")
     .def(nb::init<const std::string&>())
     .def("eval_sat", &sa::eval_sat, nb::kw_only(), "T"_a, "k"_a, "q"_a)
     .def("eval_sat_many", &sa::eval_sat_many<tensor1d>, nb::kw_only(), "T"_a, "k"_a, "q"_a, "y"_a)
     .def("get_yval", &sa::get_yval, nb::kw_only(), "T"_a, "q"_a, "k"_a)
     .def("get_yval_many", &sa::get_yval_many<tensor1d>, nb::kw_only(), "T"_a, "k"_a, "q"_a, "y"_a)
     .def("get_approx1d", &sa::get_approx1d, nb::kw_only(), "k"_a, "q"_a, nb::rv_policy::reference)
     .def_prop_ro("invlnp", &sa::get_invlnp, nb::rv_policy::reference)
     .def("solve_for_T", &sa::solve_for_T, nb::kw_only(), "propval"_a, "k"_a, "q"_a, "bits"_a = 64, "max_iter"_a=100, "boundsftol"_a=1e-13)
     .def("get_T_from_p", &sa::get_T_from_p, nb::kw_only(), "p"_a)
     .def("get_vaporquality", &sa::get_vaporquality, nb::kw_only(), "T"_a, "propval"_a, "k"_a)
     .def("add_variable", &sa::add_variable, nb::kw_only(), "k"_a, "caller"_a)
     .def("solve_for_Tq_DX", &sa::solve_for_Tq_DX)
     .def("solve_for_Tq_DX_many", &sa::solve_for_Tq_DX_many<tensor1d>)
     .def("flash", &sa::flash)
     .def("flash_many", &sa::flash_many<tensor1d>)
    ;
    
    nb::class_<FlashPhase>(m, "FlashPhase")
        .def_ro("rho_molm3", &FlashPhase::rho_molm3)
        .def_ro("qmolar", &FlashPhase::qmolar)
        .def_ro("mole_fractions", &FlashPhase::mole_fractions);
    
    nb::class_<FlashSolution>(m, "FlashSolution")
        .def_ro("T_K", &FlashSolution::T_K)
        .def_ro("Nphases", &FlashSolution::Nphases)
        .def_ro("phases", &FlashSolution::phases);
    
    nb::class_<TrhoLookup>(m, "TrhoLookup")
        .def_ro("T", &TrhoLookup::T)
        .def_ro("rho", &TrhoLookup::rho)
        .def_ro("d2", &TrhoLookup::d2);
    
    using namespace teqpflsh::properties::interfaces;
    nb::class_<HelmholtzInterface>(m, "HelmholtzInterface");
    nb::class_<teqpHelmholtzInterface, HelmholtzInterface>(m, "teqpHelmholtzInterface")
        .def(nb::init<const std::string&, const std::string&>(), nb::kw_only(), "ideal_gas"_a, "residual"_a);
    
    using mf = teqpflsh::MainFlasher;
    nb::class_<mf>(m, "MainFlasher")
        .def(nb::init<RegionedFlasher&&, const std::optional<sa>&, const std::shared_ptr<HelmholtzInterface>>(), nb::kw_only(), "regions"_a, "superancillary"_a, "helm"_a)
      .def("flash", &mf::flash, "proppair"_a, "val1"_a, "val2"_a)
      .def("flash_many", &mf::flash_many<tensor1d>, "proppair"_a, "val1"_a, "val2"_a, "T"_a, "rho"_a, "q"_a)
      .def_prop_ro("regioned_flasher", &mf::get_regioned_flasher, nb::rv_policy::reference)
    ;
    
    m.def("get_pair_from_chars", &get_pair_from_chars);
    m.def("get_property_chars", &get_property_chars);
    m.def("get_pair_log_scaling", &get_pair_log_scaling);
    
    // Exposed to make it easier to experiment with this function
    m.def("toms748_solve", [](const std::function<double(double)>& f, double a, double b, unsigned int bits, std::size_t max_iter){
        using namespace boost::math::tools;
        auto result = toms748_solve(f, a, b, eps_tolerance<double>(bits), max_iter);
        return std::make_tuple(result, max_iter);
    });
}
