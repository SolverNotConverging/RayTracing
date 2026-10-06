#include "crt/propagation.hpp"
#include "crt/visualization.hpp"

#include <matplot/matplot.h>

#include <complex>
#include <iostream>

int main()
{
    using namespace crt;
    using namespace matplot;

    const Vec3 tx{0.0, 0.0, 0.0};
    const Vec3 rx{10.0, 4.0, 2.0};

    constexpr double frequency = 3.5e9;

    const auto path =
        freeSpacePath(tx, rx, frequency);

    std::cout
        << "Length: "
        << path.length
        << " m\n";

    std::cout
        << "Delay: "
        << path.delay * 1e9
        << " ns\n";

    std::cout
        << "h = "
        << path.coefficient
        << '\n';

    std::cout
        << "|h| = "
        << std::abs(path.coefficient)
        << '\n';

    std::cout
        << "phase = "
        << std::arg(path.coefficient)
        << " rad\n";

    // -----------------------------
    // Visualization
    // -----------------------------

    figure();

    hold(on);

    plotRaySegment(tx, rx);

    plotPoint(tx);
    plotPoint(rx);

    xlabel("x [m]");
    ylabel("y [m]");
    zlabel("z [m]");

    title("Line-of-sight propagation path");

    grid(on);

    view(45, 25);

    hold(off);

    show();
}