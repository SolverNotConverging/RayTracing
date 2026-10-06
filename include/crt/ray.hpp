#pragma once

#include <Eigen/Core>

namespace crt {

    using Vec3 = Eigen::Vector3d;

    struct Ray {
        Vec3 origin;
        Vec3 direction;

        Ray(const Vec3& origin_, const Vec3& direction_)
            : origin(origin_),
              direction(direction_.normalized())
        {
        }

        [[nodiscard]]
        Vec3 at(double t) const
        {
            return origin + t * direction;
        }
    };

} // namespace crt