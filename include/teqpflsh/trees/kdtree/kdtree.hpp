#pragma once

#include "nanoflann.hpp"

namespace teqpflsh::kdtree {

template<typename T=double>
class AbstractScaler{
public:
    virtual T scale(const T&) const = 0;
    virtual T descale(const T&) const = 0;
    virtual ~AbstractScaler() = default;
};

template<typename T=double>
class NoOpScaler : public AbstractScaler<T>{
protected:
    T xmin, xmax;
public:
    NoOpScaler(){};
    NoOpScaler(NoOpScaler&&) = default;
    T scale (const T& x) const override { return x; }
    T descale (const T& xscale) const override { return xscale; }
};

template<typename T=double>
class MinMaxScaler : public AbstractScaler<T>{
protected:
    T xmin, xmax;
public:
    MinMaxScaler(T xmin, T xmax) : xmin(xmin), xmax(xmax){}
    MinMaxScaler(MinMaxScaler&&) = default;
    
    T scale (const T& x) const override {
        return (x-xmin)/(xmax-xmin);
    }
    T descale (const T& xscaled) const override {
        return (xmax - xmin)*xscaled + xmin;
    }
};

template<typename T=double>
class MinMaxLogScaler : public AbstractScaler<T>{
protected:
    T log10xmin, log10xmax;
public:
    MinMaxLogScaler(T xmin, T xmax) : log10xmin(log10(xmin)), log10xmax(log10(xmax)){}
    MinMaxLogScaler(MinMaxLogScaler&&) = default;
    
    T scale (const T& x) const override {
        return (log10(x)-log10xmin)/(log10xmax-log10xmin);
    }
    T descale (const T& xscaled) const override {
        auto log10x = (log10xmax - log10xmin)*xscaled + log10xmin;
        return pow(10.0, log10x);
    }
};

template <typename T, typename ArrType=Eigen::Map<const Eigen::ArrayXd>>
struct ViewAdapter
{
    using self_t = ViewAdapter<T, ArrType>;
    const std::shared_ptr<const ArrType> x, y;
    const std::shared_ptr<AbstractScaler<T>> xscaler, yscaler;
    
    using coord_t = T;  //!< The type of each coordinate
    
    ViewAdapter(const std::shared_ptr<const ArrType> x,
                const std::shared_ptr<const ArrType> y,
                const std::shared_ptr<AbstractScaler<T>> &xscaler = nullptr,
                const std::shared_ptr<AbstractScaler<T>> &yscaler = nullptr)
        : x(x), y(y), xscaler(xscaler), yscaler(yscaler) {
            if ( this->x->size()!= this->y->size()){ throw std::invalid_argument("x and y arguments are not the same size"); }
        };
    std::size_t get_used_bytes() const{
        return 2*x->size()*sizeof(decltype((*x)[0])) + 2*8;
    }
    
    ViewAdapter (const ViewAdapter&) = default;
    ViewAdapter& operator= (const ViewAdapter&) = default;
    ViewAdapter(ViewAdapter&& other) : x(std::move(x)), y(std::move(y)), xscaler(std::move(xscaler)), yscaler(std::move(yscaler)) {};

    // The number of data points
    inline size_t kdtree_get_point_count() const { return x->size(); }

    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0)
            return (xscaler) ? xscaler->scale((*x)(idx)) : (*x)(idx);
        else
            return (yscaler) ? yscaler->scale((*y)(idx)) : (*y)(idx);
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
        return false;
    }
};

template<typename ArrType=Eigen::ArrayXd>
class L2Tree{

public:

    using num_t = double; // TODO: investigate making single precision instead of double
    const ViewAdapter<double, ArrType> cloud;
    using pc2 = decltype(cloud);
    
    // Define the kd-tree index along with its distance metric
    using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<num_t, pc2>,
        pc2,  // The adapter class
        2 /* dim */
    >;

    my_kd_tree_t m_index;

    L2Tree(const std::shared_ptr<const ArrType> &x,
           const std::shared_ptr<const ArrType> &y,
           int tree_depth,
           const std::shared_ptr<AbstractScaler<num_t>> &xscaler = nullptr, // should be unique_ptr probably
           const std::shared_ptr<AbstractScaler<num_t>> &yscaler = nullptr):
    cloud(x, y, xscaler, yscaler), m_index(2 /* dim */, cloud, tree_depth){}
          
    L2Tree (const L2Tree&) = default;
    L2Tree (L2Tree&&) = default;
    L2Tree& operator= (const L2Tree&) = delete;
    
    std::size_t get_used_bytes(){
        return cloud.get_used_bytes() + m_index.usedMemory(m_index);
    }
          
    inline auto scale_pt(const Eigen::Array2d& query_pt) const {
        return (EArray2() <<
                  ((cloud.xscaler) ? cloud.xscaler->scale(query_pt[0]) : query_pt[0]),
                  ((cloud.yscaler) ? cloud.yscaler->scale(query_pt[1]) : query_pt[1])
                ).finished();
    }
    
    std::tuple<std::size_t, double> get_nearest_indexd2(const double x, const double y) const{
        // do a knn search
        const size_t                   num_results = 1;
        size_t                         ret_index;
        num_t                          out_dist_sqr;
        nanoflann::KNNResultSet<num_t> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr);
        Eigen::Array2d scaled_pt = (EArray2() <<
                  ((cloud.xscaler) ? cloud.xscaler->scale(x) : x),
                  ((cloud.yscaler) ? cloud.yscaler->scale(y) : y)
        ).finished();
        m_index.findNeighbors(resultSet, scaled_pt.data(), nanoflann::SearchParameters(10));
        return std::make_tuple(ret_index, out_dist_sqr);
    }

    std::tuple<std::size_t, double> get_nearest_indexd2(const Eigen::Array2d& query_pt) const{
        // do a knn search
        const size_t                   num_results = 1;
        size_t                         ret_index;
        num_t                          out_dist_sqr;
        nanoflann::KNNResultSet<num_t> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr);
        auto scaled_pt = scale_pt(query_pt);
        m_index.findNeighbors(resultSet, scaled_pt.data(), nanoflann::SearchParameters(10));
        return std::make_tuple(ret_index, out_dist_sqr);
    }
    
    template<typename Container, typename ContainerInt>
    void get_nearest_indexd2_many(const Container& x, const Container &y, ContainerInt &idx, Container &d2) const{
        for (auto i = 0; i < x.size(); ++i){
            auto [idx_, d2_] = get_nearest_indexd2(x(i), y(i));
            idx(i) = static_cast<int>(idx_);
            d2(i) = d2_;
        }
    }

//    template<typename T>
//    auto get_nearest_indexd2_T(const T& query_pt) const{
//        // do a knn search
//        const size_t                   num_results = 1;
//        size_t                         ret_index;
//        num_t                          out_dist_sqr;
//        nanoflann::KNNResultSet<num_t> resultSet(num_results);
//        resultSet.init(&ret_index, &out_dist_sqr);
//        auto scaled_pt = scale_pt(query_pt);
//        m_index.findNeighbors(resultSet, scaled_pt.data(), nanoflann::SearchParameters(10));
//        return std::make_tuple(ret_index, out_dist_sqr);
//    }

};

} // namespace teqpflsh::kdtree
