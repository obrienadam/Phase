#include "Intersections.h"

std::vector<Point3D> intersections(const Sphere &sphere, const LineSegment3D &ln)
{
    Point3D x0 = ln.ptA() - sphere.centroid();
    Vector3D r = ln.r();

    Scalar a = r.magSqr();
    Scalar b = 2. * dot(x0, r.magSqr());
    Scalar c = x0.magSqr() - std::pow(sphere.radius(), 2);

    Scalar disc = b * b - 4. * a * c;

    std::vector<Point3D> result;

    if(disc < 0.)
        return std::vector<Point3D>();
    else if (disc == 0.)
    {
        Scalar t = -b / (2. * a);

        if(t >= 0. && t <= 1.)
            result.push_back(ln.ptA() + t * r);
    }
    else
    {
        disc = std::sqrt(disc);

        Scalar t1 = (-b - disc) / (2. * a);
        Scalar t2 = (-b + disc) / (2. * a);

        result.reserve(2);

        if(t1 >= 0. && t1 <= 1.)
            result.push_back(ln.ptA() + t1 * r);

        if(t2 >= 0. && t2 <= 1.)
            result.push_back(ln.ptA() + t2 * r);
    }

    return result;
}

std::vector<Point3D> intersections(const Sphere &sphere, const Ray3D &ray)
{
    Point3D x0 = ray.x0() - sphere.centroid();

    Scalar a = ray.r().magSqr();
    Scalar b = 2. * dot(x0, ray.r());
    Scalar c = x0.magSqr() - std::pow(sphere.radius(), 2);

    Scalar disc = b * b - 4. * a * c;

    std::vector<Point3D> result;

    if(disc < 0.)
        return std::vector<Point3D>();
    else if (disc == 0.)
    {
        Scalar t = -b / (2. * a);

        if(t >= 0.)
            result.push_back(ray(t));
    }
    else
    {
        disc = std::sqrt(disc);

        Scalar t1 = (-b - disc) / (2. * a);
        Scalar t2 = (-b + disc) / (2. * a);

        result.reserve(2);

        if(t1 >= 0.)
            result.push_back(ray(t1));

        if(t2 >= 0.)
            result.push_back(ray(t2));
    }

    return result;
}

std::pair<Point3D, bool> intersection(const Plane &plane, const Ray3D &ray)
{
    Scalar t = (plane.coeffs()[3] - dot(plane.n(), ray.x0())) / dot(plane.n(), ray.r());
    return std::make_pair(ray(t), t >= 0.);
}

std::pair<Point3D, bool> intersection(const Plane &plane, const LineSegment3D &ln)
{
    Scalar t = (plane.coeffs()[3] - dot(plane.n(), ln.ptA())) / dot(plane.n(), ln.r());
    return std::make_pair(ln.ptA() + t * ln.r(), t >= 0. && t <= 1.);
}
