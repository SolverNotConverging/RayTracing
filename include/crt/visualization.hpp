#pragma once

#include "crt/ray.hpp"

#include <matplot/matplot.h>

#include <vector>

namespace crt {

    inline void plotRaySegment(
        const Vec3& start,
        const Vec3& end)
    {
        using namespace matplot;

        std::vector<double> x{
            start.x(),
            end.x()
        };

        std::vector<double> y{
            start.y(),
            end.y()
        };

        std::vector<double> z{
            start.z(),
            end.z()
        };

        plot3(x, y, z)->line_width(2.0);
    }

    inline void plotPoint(const Vec3& p)
    {
        using namespace matplot;

        scatter3(
            std::vector<double>{p.x()},
            std::vector<double>{p.y()},
            std::vector<double>{p.z()}
        );
    }

} // namespace crt