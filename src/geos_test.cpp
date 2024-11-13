/*
* # GEOS C++ example 1
*
* Reads two WKT representations and calculates the
* intersection, prints it out, and cleans up.
*
* In general, to avoid API changes, stick to operations
* on Geometry. The more esoteric APIs are more likely
* to change between versions.
*/

#include <iostream>

/* For geometry operations */
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>

/* For WKT read/write */
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>

/* Geometry/GeometryFactory */
using namespace geos::geom;

/* WKTReader/WKTWriter */
using namespace geos::io;

int main()
{
    /* New factory with default (float) precision model */
    GeometryFactory::Ptr factory = GeometryFactory::create();

    /*
    * Reader requires a factory to bind the geometry to
    * for shared resources like the PrecisionModel
    */
    WKTReader reader(*factory);

    /* Input WKT strings */
    std::string wkt_a("POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))");
    std::string wkt_b("POLYGON((5 5, 15 5, 15 15, 5 15, 5 5))");

    /* Convert WKT to Geometry */
    std::unique_ptr<Geometry> geom_a(reader.read(wkt_a));
    std::unique_ptr<Geometry> geom_b(reader.read(wkt_b));

    /* Calculate intersection */
    std::unique_ptr<Geometry> inter = geom_a->intersection(geom_b.get());

    /* Convert Geometry to WKT */
    WKTWriter writer;
    writer.setTrim(true); /* Only needed before GEOS 3.12 */
    std::string inter_wkt = writer.write(inter.get());

    /* Print out results */
    std::cout << "Geometry A:         " << wkt_a << std::endl;
    std::cout << "Geometry B:         " << wkt_b << std::endl;
    std::cout << "Intersection(A, B): " << inter_wkt << std::endl;
    
    auto create_circle = [&factory](int N, double x0, double y0, double r){
        CoordinateSequence seq;
        const double PI = 3.141592653589793238462643383279502884L;
        for (auto t = 0; t <= 2*PI; t += 2*PI/(N-1)){
            seq.add(x0+r*cos(t), y0+r*sin(t));
        }
        return std::unique_ptr<Polygon>(factory->createPolygon(std::move(seq)));
    };
    auto circle_a = create_circle(100, 0, 0, 1);
    auto circle_b = create_circle(100, 0.5, 1, 1);
    
    std::unique_ptr<Geometry> inter2 = circle_a->intersection(circle_b.get());
    
    for (auto ipoly = 0; ipoly < inter2->getNumGeometries(); ipoly++){
        inter->getGeometryN(ipoly)->clone();
    }
    std::cout << inter2->toString() << std::endl;

}
