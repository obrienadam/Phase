#include "Intersection.h"
#include "Math/StaticMatrix.h"

std::pair<Point2D, bool> intersection(const Line2D &lnA, const Line2D &lnB)
{
    auto A = StaticMatrix<2, 2>({
                                    lnA.d().x, -lnB.d().x,
                                    lnA.d().y, -lnB.d().y
                                });

    auto b = StaticMatrix<2, 1>({
                                    -lnA.r0().x + lnB.r0().x,
                                    -lnA.r0().y + lnB.r0().y
                                });

    A.solve(b);

    return std::make_pair(lnA(b(0, 0)), true);
}

std::pair<Point2D, bool> intersection(const Ray2D &ray, const LineSegment2D &ln)
{
    auto A = StaticMatrix<2, 2>({
                                    ray.r().x, -(ln.ptB().x - ln.ptA().x),
                                    ray.r().y, -(ln.ptB().y - ln.ptA().y)
                                });

    auto b = StaticMatrix<2, 1>({
                                    -ray.x0().x + ln.ptA().x,
                                    -ray.x0().y + ln.ptA().y
                                });

    A.solve(b);

    Point2D xc = ray(b(0, 0));

    return std::make_pair(xc, ray.isBounded(xc) && ln.isBounded(xc));
}

StaticVector<Point2D, 2> intersection(const Circle &circle, const Ray2D &ray)
{
    Vector2D r = ray.x0() - circle.centroid();

    Scalar a = ray.r().magSqr();
    Scalar b = 2 * dot(r, ray.r());
    Scalar c = r.magSqr() - circle.radius() * circle.radius();
    Scalar disc = b * b - 4 * a * c;

    StaticVector<Point2D, 2> result;

    if (disc == 0. && b < 0.)
    {
        result.push_back(ray(-b / (2. * a)));
    }
    else if (disc > 0.)
    {
        Scalar t1 = (-b - sqrt(disc)) / (2 * a);
        Scalar t2 = (-b + sqrt(disc)) / (2 * a);

        if (t1 > 0.)
            result.push_back(ray(t1));

        if (t2 > 0.)
            result.push_back(ray(t2));
    }

    return result;
}
