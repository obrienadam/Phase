#include "Tensor3D.h"

Tensor3D &Tensor3D::operator +=(const Tensor3D &tau)
{
    xx += tau.xx;
    xy += tau.xy;
    xz += tau.xz;
    yx += tau.yx;
    yy += tau.yy;
    yz += tau.yz;
    zx += tau.zx;
    zy += tau.zy;
    zz += tau.zz;
    return *this;
}

Tensor3D &Tensor3D::operator -=(const Tensor3D &tau)
{
    xx -= tau.xx;
    xy -= tau.xy;
    xz -= tau.xz;
    yx -= tau.yx;
    yy -= tau.yy;
    yz -= tau.yz;
    zx -= tau.zx;
    zy -= tau.zy;
    zz -= tau.zz;
    return *this;
}

Tensor3D &Tensor3D::operator *=(Scalar a)
{
    xx *= a;
    xy *= a;
    xz *= a;
    yx *= a;
    yy *= a;
    yz *= a;
    zx *= a;
    zy *= a;
    zz *= a;
    return *this;
}

Tensor3D &Tensor3D::operator /=(Scalar a)
{
    xx /= a;
    xy /= a;
    xz /= a;
    yx /= a;
    yy /= a;
    yz /= a;
    zx /= a;
    zy /= a;
    zz /= a;
    return *this;
}

Tensor3D Tensor3D::transpose() const
{
    return Tensor3D(xx, yx, zx,
                    xy, yy, zy,
                    xz, yz, zz);
}

//- Functions
Vector3D dot(const Tensor3D &tau, const Vector3D &u)
{
    return Vector3D(
                tau.xx * u.x + tau.xy * u.y + tau.xz * u.z,
                tau.yx * u.x + tau.yy * u.y + tau.yz * u.z,
                tau.zx * u.x + tau.zy * u.y + tau.zz * u.z
                );
}

Tensor3D outer(const Vector3D &u, const Vector3D &v)
{
    return Tensor3D(
                u.x * v.x, u.x * v.y, u.x * v.z,
                u.y * v.x, u.y * v.y, u.y * v.z,
                u.z * v.x, u.z * v.y, u.z * v.z
                );
}
