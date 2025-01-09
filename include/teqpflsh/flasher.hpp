#pragma once
#include <string>
#include <random>
#include <concepts>

#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/algorithms/iteration.hpp"

#include "teqpflsh/polyquadtree.hpp"
#include "teqpflsh/trees/kdtree/kdtree.hpp"
#include "teqpflsh/properties/properties.hpp"
#include "teqpflsh/properties/interfaces.hpp"
#include "teqpflsh/superancillary/superancillary.hpp"

#include "geos/triangulate/polygon/ConstrainedDelaunayTriangulator.h"
#include "geos/triangulate/polygon/PolygonTriangulator.h"
#include "geos/geom/prep/PreparedGeometryFactory.h"
#include "geos/operation/valid/MakeValid.h"
#include "geos/operation/valid/IsSimpleOp.h"

namespace teqpflsh {

using ArrayType = Eigen::ArrayXd;


/**
 Derived from: https://stackoverflow.com/a/47418580
 */
template<typename Geo, typename Generator>
auto point_in_triangle(const Geo&tri, Generator& gen){
    std::uniform_real_distribution<> uni;
    auto co = tri->getCoordinates();
    
    auto r1 = uni(gen), r2 = uni(gen);
    auto q = std::abs(r1 - r2);
    auto s = q;
    auto t = 0.5 * (r1 + r2 - q);
    auto u = 1 - 0.5 * (q + r1 + r2);
    auto x = s * co->getX(0) + t * co->getX(1) + u * co->getX(2);
    auto y = s * co->getY(0) + t * co->getY(1) + u * co->getY(2);
    return std::make_tuple(x, y);
}

template <class T, class = void>
struct has_operatorbracket1 : std::false_type {};

template <class T>
struct has_operatorbracket1 < T,
            std::void_t<decltype(std::declval<T>()[0])>> :
            std::true_type {};

/// Sample a given number of points that are within the polygon by first doing triangulation of
/// the bounding polygon with triangulation, sampling the triangles according to their relative
/// area, and then randomly sampling a point inside a triangle
template<typename Container>
void sample_random(const Geometry* geo, std::size_t Nsamples, Container& destx, Container& desty){
    if(destx.size() != Nsamples){ throw std::invalid_argument("destx size does not equal Nsamples"); }
    if(desty.size() != Nsamples){ throw std::invalid_argument("desty size does not equal Nsamples"); }
    
    auto tri = geos::triangulate::polygon::PolygonTriangulator::triangulate(geo);
    
    // Get the areas for each triangle
    auto Ntri = tri->getNumGeometries();
    std::vector<double> areas(Ntri);
    for (auto igeo = 0U; igeo < Ntri; ++igeo){
        areas[igeo] = tri->getGeometryN(igeo)->getArea();
    }
    // Prepare the weighting function
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(areas.begin(), areas.end());
    
    for (auto rngcounter = 0U; rngcounter < Nsamples; ++rngcounter){
        // Randomly select a triangle within the triangulation weighted according to
        // its relative area.
        auto selected_triangle = tri->getGeometryN(d(gen));
        // And select a point within this triangle
        auto [x, y] = point_in_triangle(selected_triangle, gen);
        if constexpr(has_operatorbracket1<Container>::value){
            destx[rngcounter] = x;
            desty[rngcounter] = y;
        }
        else{
            destx(rngcounter) = x;
            desty(rngcounter) = y;
        }
    }
}

/**
 This class is a generic type that represents a 2D region defined by its bounding envelope
 which is expected to be a closed & non-intersecting polygon. A quadtree can be constructed covering the
 entire region so that inside/outside determination is very fast (order of nanoseconds).  The polygon is segmented
 into parts that are mostly rectangles coincident with the rectangular grid of the quadtree.
 */
class QuadRegion2D{
public:
    struct Envelope{
        double x_min, y_min, x_max, y_max;
    };
    
private:
    auto build_envelope(const std::unique_ptr<Geometry>& poly){
        if (poly == nullptr){
            throw std::invalid_argument("Bad polygon");
        }
        auto envelope = poly->getEnvelope();
        const auto coords = envelope->getCoordinates();
        
        // 5 points, starting in lower left and going counter-clockwise
        // so max is at index 2
        auto x_min = coords->getX(0);
        auto y_min = coords->getY(0);
        auto x_max = coords->getX(2);
        auto y_max = coords->getY(2);
        
        return Envelope{x_min, y_min, x_max, y_max};
    }
    template<typename Vector>
    auto create_base_polygon(const Vector&x, const Vector&y){
        CoordinateSequence seq;
        for (auto i=0; i < x.size(); ++i){
            seq.add(x[i], y[i]);
        }
        seq.closeRing();
        auto ptr = std::unique_ptr<Geometry>(factory->createPolygon(std::move(seq)));
        return ptr;
    }
    GeometryFactory::Ptr factory;
    std::unique_ptr<Geometry> m_bounding_polygon;
    const Envelope envelope;
    PolyQuadTree quadtree;
    std::unique_ptr<geos::geom::prep::PreparedGeometry> m_preppoly;
    
private:
    QuadRegion2D(std::unique_ptr<Geometry>&& bounding_polygon) : m_bounding_polygon(std::move(bounding_polygon)), envelope(build_envelope(m_bounding_polygon)), quadtree(m_bounding_polygon->getCoordinates()), m_preppoly(geos::geom::prep::PreparedGeometryFactory().create(m_bounding_polygon.get()))
    {}
public:
    
    /// Public constructor taking the x and y coordinates of the geometry
    QuadRegion2D(const std::vector<double>&x, const std::vector<double>& y) : factory(GeometryFactory::create()), m_bounding_polygon(create_base_polygon(x, y)), envelope(build_envelope(m_bounding_polygon)), quadtree(m_bounding_polygon->getCoordinates()), m_preppoly(geos::geom::prep::PreparedGeometryFactory().create(m_bounding_polygon.get()))
    {
        // TODO: see IsSimpleOp to identify issues in polygon at construction time
    }
    /// Public constructor taking the x and y coordinates of the geometry as Eigen::Array
    QuadRegion2D(const Eigen::ArrayXd&x, const Eigen::ArrayXd& y) : factory(GeometryFactory::create()), m_bounding_polygon(create_base_polygon(x, y)), envelope(build_envelope(m_bounding_polygon)), quadtree(m_bounding_polygon->getCoordinates()), m_preppoly(geos::geom::prep::PreparedGeometryFactory().create(m_bounding_polygon.get())) {
    }
    QuadRegion2D(const QuadRegion2D&) = delete;
    auto get_envelope() const {
        return envelope;
    }
    bool in_envelope(double x, double y) const {
        auto& e = envelope;
        return (x >= e.x_min && x <= e.x_max && y >= e.y_min && y <= e.y_max);
    }
    bool contains(double x, double y) const {
        auto pt = factory->createPoint(CoordinateXY{x, y});
        return m_preppoly->contains(pt.get());
    }
    bool distance(double x, double y) const {
        auto pt = factory->createPoint(CoordinateXY{x, y});
        return m_preppoly->distance(pt.get());
    }
        
    auto get_coords_xy(){
        std::vector<CoordinateXY> coordsxy;
        m_bounding_polygon->getCoordinates()->toVector(coordsxy);
        std::vector<double> x, y;
        for (auto&&xy : coordsxy){
            x.push_back(xy.x);
            y.push_back(xy.y);
        }
        return std::make_tuple(x, y);
    }
    auto get_nonsimple_xy(){
        geos::operation::valid::IsSimpleOp op(m_bounding_polygon.get());
        op.setFindAllLocations(true);
        std::vector<CoordinateXY> coordsxy = op.getNonSimpleLocations();
        std::vector<double> x, y;
        for (auto&&xy : coordsxy){
            x.push_back(xy.x);
            y.push_back(xy.y);
        }
        return std::make_tuple(x, y);
    }
    Geometry* get_bounding_polygon(){
        return m_bounding_polygon.get();
    }
    
    void do_splits(int Nsplits /*, constraints*/){
        quadtree.do_splits(Nsplits/*, constraints*/);
    }
    const PolyQuadTree& get_quadtree_ro() const{
        return quadtree;
    }
    PolyQuadTree& get_quadtree_rw(){
        return quadtree;
    }
    
    /// Do Delaunay triangulation which tends to yield well-shaped triangles
    auto do_Delaunay_triangulation(){
        return geos::triangulate::polygon::ConstrainedDelaunayTriangulator::triangulate(m_bounding_polygon.get());
    }
    /// Do "fast" triangulation which tends to yield less-well-shaped triangles but is faster than Delaunay triangulation
    auto do_fast_triangulation(){
        return geos::triangulate::polygon::PolygonTriangulator::triangulate(m_bounding_polygon.get());
    }
    
    /**
     Derived from: https://stackoverflow.com/a/47418580
     */
    template<typename Geo, typename Generator>
    auto point_in_triangle(const Geo&tri, Generator& gen){
        std::uniform_real_distribution<> uni;
        auto co = tri->getCoordinates();
        
        auto r1 = uni(gen), r2 = uni(gen);
        auto q = std::abs(r1 - r2);
        auto s = q;
        auto t = 0.5 * (r1 + r2 - q);
        auto u = 1 - 0.5 * (q + r1 + r2);
        auto x = s * co->getX(0) + t * co->getX(1) + u * co->getX(2);
        auto y = s * co->getY(0) + t * co->getY(1) + u * co->getY(2);
        return std::make_tuple(x, y);
    }
    
    template <class T, class = void>
    struct has_operatorbracket1 : std::false_type {};

    template <class T>
    struct has_operatorbracket1 < T,
                std::void_t<decltype(std::declval<T>()[0])>> :
                std::true_type {};
    
    /// Sample a given number of points that are within the polygon by first doing triangulation of
    /// the bounding polygon with triangulation, sampling the triangles according to their relative
    /// area, and then randomly sampling a point inside a triangle
    template<typename Container>
    void sample_random(std::size_t Nsamples, Container& destx, Container& desty){
        return teqpflsh::sample_random<Container>(m_bounding_polygon.get(), Nsamples, destx, desty);
    }
    
    /**
    Obtain a tensor grid of points that are within the bounding polygon of the region
     
     The number of points that are within the bounding polygon is not known until after the function completes
     so to be sure, allocate buffers as large as the product of the dimensions of the gridded points
     */
    template<typename Container>
    std::size_t sample_gridded(const Container& gridx, const Container& gridy, Container& bufx, Container& bufy){
        
        if (bufx.size() != bufy.size()){ throw std::invalid_argument("Output buffers must be the same size"); }
        
        std::size_t index = 0;
        for (auto i = 0U; i < gridx.size(); ++i){
            for (auto j = 0U; j < gridy.size(); ++j){
                auto pt = factory->createPoint(CoordinateXY{gridx(i), gridy(j)});
                if (m_preppoly->contains(pt.get())){
                    bufx(index) = gridx(i);
                    bufy(index) = gridy(j);
                    index++;
                    if (index == bufx.size()){
                        throw std::invalid_argument("output buffer is exhausted, allocate a larger buffer");
                    }
                }
            }
        }
        return index;
    }
    
    /**
    Obtain a tensor grid of points that are within the bounding polygon of the region
     
     The number of points that are within the bounding polygon is not known until after the function completes
     so to be sure, allocate buffers as large as the product of the dimensions of the gridded points
     */
    template<typename Container>
    std::size_t sample_gridded_w_tree(const Container& gridx, const Container& gridy, Container& bufx, Container& bufy){
        if (bufx.size() != bufy.size()){ throw std::invalid_argument("Output buffers must be the same size"); }
        auto& tree = get_quadtree_rw();
        
        std::size_t index = 0;
        for (auto i = 0U; i < gridx.size(); ++i){
            for (auto j = 0U; j < gridy.size(); ++j){
                auto& node = tree.getNode(gridx(i), gridy(j), false);
                auto status = node.get_contents().status;
                
                if (status == PQTStatus::outside){
                    // Node is completely outside this sample is ignored
                    continue;
                }
                else if (status==PQTStatus::inside){
                    // Node is completely inside and kept
                    bufx(index) = gridx(i);
                    bufy(index) = gridy(j);
                    index++;
                    if (index == bufx.size()){
                        throw std::invalid_argument("output buffer is exhausted, allocate a larger buffer");
                    }
                }
                else{
                    // We need to use the polygon of the intersection node to determine what to do
                    auto pt = factory->createPoint(CoordinateXY{gridx(i), gridy(j)});
                    if (m_preppoly->contains(pt.get())){
                        bufx(index) = gridx(i);
                        bufy(index) = gridy(j);
                        index++;
                        if (index == bufx.size()){
                            throw std::invalid_argument("output buffer is exhausted, allocate a larger buffer");
                        }
                    }
                }
            }
        }
        return index;
    }
};

std::ostream& operator<<(std::ostream& os, const QuadRegion2D::Envelope& e)
{
    os << e.x_min << "," << e.x_max << "," << e.y_min << "," << e.y_max;
    return os;
}

//static_assert(std::is_move_constructible_v<QuadRegion2D>);
//static_assert(std::is_move_assignable_v<QuadRegion2D>);

template<typename T>
concept CallableMin = requires(T t) {
    { t.min() };
};

template<typename T>
auto minmax(const T& v) requires CallableMin<T>{
    return std::make_pair(v.min(), v.max());
}

template<typename Vector>
auto minmax(const Vector& v){
    if (v.size() == 0){
        throw std::invalid_argument("can't get min of empty");
    }
    else{
        using T = std::decay_t<decltype(v[0])>;
        T minval = v[0], maxval = v[0];
        for (auto i = 1L; i < v.size(); ++i){
            if(v[i] < minval){ minval = v[i]; }
            if(v[i] > maxval){ maxval = v[i]; }
        }
        return std::make_pair(minval, maxval);
    }
}

template<typename Container>
auto geomspace(const Container& v, Eigen::Index N){
    auto [min_, max_] = minmax(v);
    return (Eigen::ArrayXd::LinSpaced(N, log(min_), log(max_)).exp()).eval();
}

struct TrhoLookup{
    double T, rho, d2;
};

/**
 A class that contains all the transformed polygons as well as K-D trees for each input pair
 that indexes back to the T, rho coordinates used to define the K-D tree
 */
template<typename ArrayType>
struct ThermodynamicRegion{
    using QuadRegionPtr = std::shared_ptr<QuadRegion2D>;
    using L2Tree_t = kdtree::L2Tree<ArrayType>;
    using L2TreePtr = std::shared_ptr<L2Tree_t>;
    using PropertySet = teqpflsh::properties::PropertySet<ArrayType>;
    using PropertyPairs = teqpflsh::properties::PropertyPairs;
private:
    std::unordered_map<PropertyPairs, QuadRegionPtr> transformed_regions;
    std::unordered_map<PropertyPairs, L2TreePtr> kdtrees;

    PropertySet make_propset_grid(const ArrayType& z, const auto& alphamodel, Eigen::Index NT, Eigen::Index Nrho){
        using namespace std::chrono;
        auto start = high_resolution_clock::now();
        auto qr = QuadRegion2D(T_bounding, rho_bounding);
        Eigen::ArrayXd gridT = geomspace(T_bounding, NT),
                       gridrho = geomspace(rho_bounding, Nrho),
                       bufT = Eigen::ArrayXd(gridT.size()*gridrho.size()),
                       bufrho = Eigen::ArrayXd(gridrho.size()*gridT.size());
        auto Npts = qr.sample_gridded(gridT, gridrho, bufT, bufrho);
        bufT.conservativeResize(Npts);
        bufrho.conservativeResize(Npts);
        
        // And append the bounding polygon values to the points
        // to be considered when building the K-D tree
        bufT.conservativeResize(bufT.size() + T_bounding.size());
        bufrho.conservativeResize(bufrho.size() + rho_bounding.size());
        bufT.tail(T_bounding.size()) = T_bounding;
        bufrho.tail(rho_bounding.size()) = rho_bounding;
        
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
//        std::cout << duration.count() << "μs for poly stuff w/ " << Npts << std::endl;
        return {bufT, bufrho, z, alphamodel};
    }
public:
    
    const ArrayType T_bounding, ///< Bounding polygon's temperature values
                    rho_bounding; ///< Bounding polygon's density values
    const PropertySet propset_bounding;
    const PropertySet propset_Trhogrid;
    
    ThermodynamicRegion(const ArrayType &T, const ArrayType &rho, const ArrayType& z, const teqp::iteration::AlphaModel& alphamodel, Eigen::Index NT, Eigen::Index Nrho)
    : T_bounding(T), rho_bounding(rho), propset_bounding(T, rho, z, alphamodel), propset_Trhogrid(make_propset_grid(z, alphamodel, NT, Nrho)) {}
    
    // You need to worry about the rule-of-three if you use unique_ptr as the value type in map
    ThermodynamicRegion(const ThermodynamicRegion&) = default;
    ThermodynamicRegion(ThermodynamicRegion&&) = default;
    
    /// The method to actually add a L2Tree to the map
    auto add_L2Tree(PropertyPairs pair){
        using namespace kdtree;
        const auto use_log = get_pair_log_scaling(pair);
        
        auto arrays = propset_Trhogrid.get_array_ptrs(pair);
        const auto& ptsx = std::get<0>(arrays);
        const auto& ptsy = std::get<1>(arrays);
        std::shared_ptr<AbstractScaler<double>> xscaler, yscaler;
        auto [xmin, xmax] = minmax(*ptsx);
        auto [ymin, ymax] = minmax(*ptsy);
        if (std::get<0>(use_log)){
            xscaler = std::make_shared<MinMaxLogScaler<double>>(xmin, xmax);
        }
        else{
            xscaler = std::make_shared<MinMaxScaler<double>>(xmin, xmax);
        }
        if (std::get<1>(use_log)){
            yscaler = std::make_shared<MinMaxLogScaler<double>>(ymin, ymax);
        }
        else{
            yscaler = std::make_shared<MinMaxScaler<double>>(ymin, ymax);
        }
        L2TreePtr treeptr = std::make_shared<L2Tree_t>(ptsx, ptsy, 10, xscaler, yscaler);
        kdtrees.insert({pair, std::move(treeptr)});
    }
    auto remove_L2Tree(PropertyPairs pair){
        kdtrees.erase(pair);
    }
    void add_pair(const PropertyPairs pair, int Nsplit, bool and_kdtree=true){
        auto arrays = propset_bounding.get_arrays(pair);
        transformed_regions.emplace(pair, std::make_shared<QuadRegion2D>(std::get<0>(arrays), std::get<1>(arrays)));
        if (and_kdtree){
            add_L2Tree(pair);
        }
        transformed_regions.at(pair)->do_splits(Nsplit);
    }

    auto has_pair(PropertyPairs pair) const { return transformed_regions.contains(pair); }
    
    std::optional<TrhoLookup> get_starting_Trho(PropertyPairs pair, double val1, double val2) const{
        
        const auto& kd = kdtrees.at(pair);
        auto [idx, d2] = kd->get_nearest_indexd2(val1, val2);
        auto [T, rho] = propset_Trhogrid.get_Trho(idx);
        return TrhoLookup{T, rho, d2};
//
//        const auto& tr = transformed_regions.at(pair);
//        // First check whether even in envelope (bounding box) of the transformed polygon
//        if (!tr->in_envelope(val1, val2)){
//            // If not, no need to walk the tree to look for node, no solutions found
//            return std::nullopt;
//        }
//        else{
//            // Now we walk the tree to see if we are a) inside b) outside c) intersection
//            auto & pqt = tr->get_quadtree_rw();
//            const auto& node = pqt.getNode(val1, val2);
//            if (pqt.is_complete(node) || pqt.is_intersection(node)){
//                
//            }
//            else{
//                return std::nullopt;
//            }
//        }
    }
    
    
    template<typename Container>
    void get_starting_Trho_many(PropertyPairs pair, const Container& val1, const Container& val2, Container& T, Container& rho, Container& d2) const {
        for (auto i = 0; i < val1.size(); ++i){
            try{
                auto o = get_starting_Trho(pair, val1(i), val2(i));
                if (o){
                    auto& v = o.value();
                    T(i) = v.T;
                    rho(i) = v.rho;
                    d2(i) = v.d2;
                }
                else{
                    T(i) = -2;
                    rho(i) = -2;
                    d2(i) = -2;
                }
            }catch(std::exception&e){
                std::cout << e.what() << std::endl;
                T(i) = -1;
                rho(i) = -1;
                d2(i) = -1;
            }
        }
    }
    
    auto& get_propset_bounding() const { return propset_bounding; }
    auto& get_propset_Trhogrid() const { return propset_Trhogrid; }
    auto& get_transformed_regions() const { return transformed_regions; }
    const auto& get_transformed_region(PropertyPairs pair) const { return transformed_regions.at(pair); }
    const auto& get_kdtree(PropertyPairs pair) const { return kdtrees.at(pair); }
};
//static_assert(std::is_move_insertable_v<ThermodynamicRegion<Eigen::ArrayXd>>);
static_assert(std::is_move_constructible_v<ThermodynamicRegion<Eigen::ArrayXd>>);
//static_assert(std::is_move_constructible_v<ThermodynamicRegion<std::vector<double>>>);
//static_assert(std::is_move_assignable_v<ThermodynamicRegion<std::vector<double>>>);

using namespace teqp::cppinterface;

struct RegionedFlashReturn{
    double T, rho;
    int step_count;
    double maxabsr;
    teqp::iteration::StoppingConditionReason reason;
    std::string msg;
    double newton_duration_us, total_duration_us, candidate_duration_us;
};

class RegionedFlasher{
//    using ArrayType = std::vector<double>;
    using ArrayType = Eigen::ArrayXd;
    using ThermodynamicRegion_t = ThermodynamicRegion<ArrayType>;
private:
    const teqp::iteration::AlphaModel alphamodel;
    const ArrayType mole_fractions;
    const std::vector<std::shared_ptr<teqp::iteration::StoppingCondition>> stopping_conditions;
    std::vector<ThermodynamicRegion_t> regions;
    
    /// Check if a pair is ready to be used to do a flash calculation
    bool pair_is_available(properties::PropertyPairs pair) const {
        if (regions.size() == 0){
            throw std::invalid_argument("At least one region is required");
        }
        return regions[0].has_pair(pair);
    }
    std::size_t path_integration_steps = 10;
public:
    
    RegionedFlasher(const std::string&ideal_gas, const std::string& resid, const ArrayType&mole_fractions) :
        alphamodel(teqp::iteration::AlphaModel{make_model(nlohmann::json::parse(ideal_gas)), make_model(nlohmann::json::parse(resid))}),
    mole_fractions(mole_fractions),
    stopping_conditions({
        std::make_shared<teqp::iteration::MaxAbsErrorCondition>(1e-8),
        std::make_shared<teqp::iteration::MinRelStepsizeCondition>(1e-15),
        std::make_shared<teqp::iteration::NanXDXErrorCondition>(),
        std::make_shared<teqp::iteration::NegativeXErrorCondition>(),
    })
    {}
    RegionedFlasher(const RegionedFlasher&) = delete;
    RegionedFlasher(RegionedFlasher&& other)
        : alphamodel(std::move(other.alphamodel)),
          mole_fractions(std::move(other.mole_fractions)),
          stopping_conditions(std::move(other.stopping_conditions)),
          regions(std::move(other.regions)){};
    
    /// Get the mole fractions in use
    auto get_mole_fractions() const { return mole_fractions; }
    
    /** 
     Add a region in which the \f$T,\rho\f$ coordinates of the outer edge are defined
     
     \param T The array of temperatures, in K
     \param rho The array of densities, in mol/m3
     \param NT The number of points in temperature direction to generate
     \param Nrho The number of points in the density direction to generate
     */
    void add_region(const ArrayType &T, const ArrayType &rho, std::size_t NT, std::size_t Nrho){
        ThermodynamicRegion_t reg(T, rho, mole_fractions, alphamodel, NT, Nrho);
        regions.push_back(std::move(reg));
    }
    
    /**
     Remove all regions
     */
    void remove_all_regions(){
        regions.clear();
    }
    
    /// Get a const view of the regions allowing for read-only access to the regions and their subobjects
    const auto& get_regions_ro() const {
        return regions;
    }
    /// Get a writeable view of the regions allowing for read-write access to the regions and their subobjects
    auto& get_regions_rw() {
        return regions;
    }
    auto get_quadtree_intersections(properties::PropertyPairs key, double x, double y){
        std::vector<QuadRegion2D*> intersections;
        for (auto& reg: regions){
            for (auto& [k, val]: reg.get_transformed_regions()){
                if (key == k){
                    auto& qt = val->get_quadtree_rw();
                    // Check if in the envelope (bounding box) of the region
                    // This is a very fast check (two double comparisons) and prunes regions that are not candidates
                    if (qt.in_envelope(x, y)){
                        // Check if a node can be found at this location (could still be inside or outside the polygon)
                        auto& node = qt.getNode(x, y);
                        if (node.get_contents().poly){ // So it is either bisected by the polygon or entirely inside the polygon
                            intersections.push_back(val.get());
                        }
                    }
                }
            }
        }
        return intersections;
    }
    
    /// Get starting values for temperature and density
    ///
    auto get_starting_Trho(properties::PropertyPairs key, double val1, double val2) const{
        std::vector<std::tuple<const ThermodynamicRegion_t*, TrhoLookup>> intersections;
        for (const auto& reg: regions){
            auto optTrho = reg.get_starting_Trho(key, val1, val2);
            if (optTrho){
                intersections.emplace_back(&reg, std::move(optTrho.value()));
            }
        }
        return intersections;
    }
    auto get_NRIterator(const std::vector<char>& vars, const EArray2& target, double T, double rho, const Eigen::Ref<Eigen::ArrayXd> &rz, const std::tuple<bool, bool> &relative_errors, const std::vector<std::shared_ptr<teqp::iteration::StoppingCondition>>& stopping_conditions){

        return teqp::iteration::NRIterator(alphamodel, vars, target, T, rho, rz, relative_errors, stopping_conditions);
    }
    
    auto flash(const properties::PropertyPairs pair, double val1, double val2) const{
        using namespace std::chrono;
        
        if (!pair_is_available(pair)){
            throw std::invalid_argument("This pair is not available");
        }
        
        auto start_outer = high_resolution_clock::now();
        // Lookup the starting values of temperature and density;
        // and sort in terms of d^2 so the first entry will be the starting point in the most
        // likely good candidate region
        auto candidates = get_starting_Trho(pair, val1, val2);
        if (candidates.size() > 1){
            std::sort(candidates.begin(), candidates.end(), [](auto& t1, auto& t2){ return std::get<1>(t1).d2 < std::get<1>(t2).d2; });
            //        if (candidates.size() != 1){
            //            throw std::invalid_argument("Incorrect region count ("+std::to_string(candidates.size())+") for inputs: "+std::to_string(val1) + "," + std::to_string(val2));
            //        }
        }
        auto stop_cand = high_resolution_clock::now();
        auto candidate_duration_us = duration_cast<nanoseconds>(stop_cand - start_outer).count()/1000.0;
        
        for (auto& [thermregion, cand] : candidates){
            
            // Do the Newton iteration (and time it individually)
            auto start = high_resolution_clock::now();
            // Make the iterator (this copies quite a few things into the iterator)
            auto iter = teqp::iteration::NRIterator(alphamodel, get_vars(pair), (EArray2() << val1, val2).finished(), cand.T, cand.rho, mole_fractions, get_relative_errors(pair), stopping_conditions);
            iter.verbose = false;
            auto reason = iter.take_steps(20);
            auto stop = high_resolution_clock::now();
            auto newton_duration_us = duration_cast<nanoseconds>(stop - start).count()/1000.0;
            
            // A lambda function to DRY the returning of outputs
            auto result_returner = [&](){
                RegionedFlashReturn r;
                r.T = iter.get_T();
                r.rho = iter.get_rho();
                r.step_count = iter.get_step_count();
                r.maxabsr = iter.get_maxabsr();
                r.reason = reason;
                r.msg = iter.get_message();
                r.newton_duration_us = newton_duration_us;
                r.candidate_duration_us = candidate_duration_us;
                auto stop_outer = high_resolution_clock::now();
                r.total_duration_us = duration_cast<nanoseconds>(stop_outer - start_outer).count()/1000.0;
                return r;
            };
            
            // Did we get to a good stopping condition? If so, finished
            using scr = teqp::iteration::StoppingConditionReason;
            if (reason == scr::success){
                return result_returner();
            }
            else if (path_integration_steps > 0){
                // Try to approach the solution with a few steps of path integration
                auto [T, rho, newx, newy] = iter.path_integration(cand.T, cand.rho, path_integration_steps);
                
                // Do the Newton iteration (again)
                iter.reset(T, rho);
                iter.verbose = false;
                auto reason = iter.take_steps(20);
                if (reason != scr::success){
                    throw std::invalid_argument("After path integration and Newton steps a valid solution was still not obtained");
                }
                else{
                    return result_returner();
                }
            }
            else{
                throw std::invalid_argument("Newton steps were not successful and path integration is not permitted");
            }
        }
        throw std::invalid_argument("Not possible to get here, make the compiler happy");
    }
    
    /// A vectorized version, mostly for profiling purposes in Python
    template<typename Container>
    void flash_many(properties::PropertyPairs pair, const Container& val1, const Container& val2, Container&T, Container&rho, Container&steps, Container&maxabsr, Container&newtontime, Container&candtime) const {
        for (auto i = 0; i < val1.size(); ++i){
            try{
                auto r = flash(pair, val1(i), val2(i));
                T(i) = r.T;
                rho(i) = r.rho;
                steps(i) = r.step_count;
                maxabsr(i) = r.maxabsr;
                newtontime(i) = r.newton_duration_us;
                candtime(i) = r.candidate_duration_us;
            }catch(std::exception&e){
//                std::cout << e.what() << std::endl;
                T(i) = -1;
                rho(i) = -1;
                steps(i) = -1;
                maxabsr(i) = -1;
                newtontime(i) = -1;
                candtime(i) = -1;
            }
        }
    }
};

struct FlashPhase{
    double rho_molm3; ///< Density in mol/m^3
    double qmolar = -1; ///< Vapor quality (moles in this phase relative to total number of moles)
    Eigen::ArrayXd mole_fractions; ///< Mole fractions of the components in the phase
};

struct FlashSolution{
    double T_K; ///< Temperature in K
    int Nphases; ///< Number of phases present
    std::vector<FlashPhase> phases;
    double rhobulk_molm3; ///< Bulk (overall) density in mol/m^3
    
    void finalize(){
        auto predicate = [](auto& p1, auto& p2){ return p1.rho_molm3 < p2.rho_molm3; };
        std::sort(phases.begin(), phases.end(), predicate);
        Nphases = static_cast<short>(phases.size());
    }
};

class MainFlasher{
public:
    const RegionedFlasher regioned_flasher;
    const std::optional<superancillary::SuperAncillary<Eigen::ArrayXd>> superanc;
    const std::shared_ptr<properties::interfaces::HelmholtzInterface> helm;
    
    MainFlasher(RegionedFlasher&& regions,
                const std::optional<superancillary::SuperAncillary<Eigen::ArrayXd>>& ancillary_data,
                const std::shared_ptr<properties::interfaces::HelmholtzInterface>& helm
                ) : regioned_flasher(std::move(regions)), superanc(ancillary_data), helm(helm) {}
    
    const auto& get_regioned_flasher() const { return regioned_flasher; }
    
    std::optional<FlashSolution> flash(properties::PropertyPairs proppair, double val1, double val2){
        // If it has a superancillary, use it to determine if the state
        // point is possibly two-phase. If yes, this is the end of the
        // flash routine and the calculation terminates
        if (superanc){
            auto soln = superanc.value().flash(proppair, val1, val2);
            // TODO: what if you have a superancillary for a mixture?
            if (soln){
                FlashSolution flsh;

                // All phases at the same temperature
                flsh.T_K = soln.value().T;
                
                // Use the vapor quality (moles of vapor per total moles) from the superancillary solution
                double qvap = soln.value().q;
                double rhov = superanc.value().get_yval(flsh.T_K, 1, 'D');
                double rhol = superanc.value().get_yval(flsh.T_K, 0, 'D');
                flsh.phases.push_back(FlashPhase{rhov, qvap, (Eigen::ArrayXd(1) << 1.0).finished()});
                flsh.phases.push_back(FlashPhase{rhol, 1-qvap, (Eigen::ArrayXd(1) << 1.0).finished()});
                flsh.finalize();
                return flsh;
            }
        }
        // Do the homogeneous region flashing
        auto r = regioned_flasher.flash(proppair, val1, val2);
        FlashSolution flsh;
        flsh.T_K = r.T;
        flsh.phases.push_back(FlashPhase{r.rho, 999, regioned_flasher.get_mole_fractions()});
        flsh.finalize();
        return flsh;
    }
    
    template<typename Container>
    void flash_many(properties::PropertyPairs proppair, const Container& val1, const Container& val2, Container& T, Container& rho, Container& q){
        if (std::set<std::size_t>({val1.size(), val2.size(), T.size(), rho.size(), q.size()}).size() != 1){
            throw std::invalid_argument("val1, val2, T, rho, q are not all the same size");
        }
        for (auto &reg : regioned_flasher.get_regions_ro()){
            if (!reg.has_pair(proppair)){
                throw std::invalid_argument("The desired pair is not available in at least one region");
            }
        }
        
        for (auto i = 0; i < val1.size(); ++i){
            try{
                auto optsoln = flash(proppair, val1(i), val2(i));
                if (optsoln){
                    auto& r = optsoln.value();
//                    runtime_check_that(r.phases.size() > 0);
                    T(i) = r.T_K;
//                    rho(i) = r.rhobulk_molm3;
                    q(i) = r.phases.front().qmolar;
                }
            }catch(std::exception&e){
//                std::cout << e.what() << std::endl;
                T(i) = -1;
                rho(i) = -1;
                q(i) = -1;
            }
        }
    }
};


} /* namespace teqpflsh */
