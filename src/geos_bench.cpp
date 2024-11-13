#include <iostream>

/* For geometry operations */
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/intersection/Rectangle.h>
#include <geos/operation/intersection/RectangleIntersection.h>

/* Geometry/GeometryFactory */
using namespace geos::geom;

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>

TEST_CASE("circlemath", "[bench]"){
    
    GeometryFactory::Ptr factory = GeometryFactory::create();
    
    auto create_circle = [&factory](int N, double x0, double y0, double r){
        CoordinateSequence seq;
        const double MYPI = 3.141592653589793238462643383279502884L;
        for (auto t = 0.0; t < 2*MYPI; t += 2*MYPI/(N-1)){
            seq.add(x0+r*cos(t), y0+r*sin(t));
        }
        seq.closeRing();
        return std::unique_ptr<Geometry>(factory->createPolygon(std::move(seq)));
    };
    auto create_box = [&factory](double xmid, double halfdx, double ymid, double halfdy){
        CoordinateSequence seq;
        seq.add(xmid-halfdx, ymid-halfdy);
        seq.add(xmid+halfdx, ymid-halfdy);
        seq.add(xmid+halfdx, ymid+halfdy);
        seq.add(xmid-halfdx, ymid+halfdy);
        seq.closeRing();
        return std::unique_ptr<Geometry>(factory->createPolygon(std::move(seq)));
    };
    auto create_rect = [&factory](double xmid, double halfdx, double ymid, double halfdy){
        return geos::operation::intersection::Rectangle(xmid-halfdx, ymid-halfdy, xmid+halfdx, ymid+halfdy);
    };
    auto circle_a = create_circle(1000, 0, 0, 1);
    auto circle_b = create_circle(1000, 0.5, 1, 1);
    auto box_a = create_box(1, 0.5, 1, 0.5);
    auto box_b = create_box(0.5, 0.5, 0.5, 0.5);
    
    auto rect_a = create_rect(1, 0.5, 1, 0.5);
    auto rect_b = create_rect(0.5, 0.5, 0.5, 0.5);
    
    auto A1 = circle_a->intersection(box_b.get())->getArea();
    auto A2 = geos::operation::intersection::RectangleIntersection::clip(*circle_a.get(), rect_b)->getArea();
    CHECK(A1 == A2);

//    BENCHMARK("creating circle"){
//        return create_circle(100, 0, 0, 1);
//    };
//    BENCHMARK("intersection of two circles"){
//        return circle_a->intersection(circle_b.get());
//    };
//    BENCHMARK("intersection of circle with box_a"){
//        return circle_a->intersection(box_a.get());
//    };
    BENCHMARK("intersection of circle with box_b"){
        return circle_a->intersection(box_b.get());
    };
    BENCHMARK("intersection of box with box"){
        return box_b->intersection(box_a.get());
    };
    BENCHMARK("clip of circle with rect"){
        return geos::operation::intersection::RectangleIntersection::clip(*circle_a.get(), rect_b);
    };
};
