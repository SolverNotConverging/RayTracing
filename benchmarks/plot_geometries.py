"""Render benchmark geometry schematics with PyVista/VTK (no solver run)."""
import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


def render(output):
    output.mkdir(parents=True, exist_ok=True)
    plotter = pv.Plotter(shape=(2, 3), window_size=(2100, 1200), off_screen=True)
    plotter.set_background('white', all_renderers=True)
    blue, orange, green, grey = '#2067b2', '#dc762c', '#18866b', '#748595'

    def line(points, color=blue, width=4):
        points = np.asarray(points, dtype=float)
        if points.shape[1] == 2:
            points = np.column_stack((points, np.zeros(len(points))))
        plotter.add_mesh(pv.lines_from_points(points), color=color, line_width=width,
                         lighting=False)

    def label(point, text, color='#243647'):
        plotter.add_point_labels([(*point, 0)], [text], font_size=19,
                                text_color=color, show_points=False, shape=None,
                                always_visible=True)

    def dot(point, color=blue):
        plotter.add_mesh(pv.Sphere(radius=.055, center=(*point, 0)), color=color)

    def panel(index, title, caption):
        plotter.subplot(index // 3, index % 3)
        plotter.add_text(title, position=(24, 548), font_size=18, color='#172e44')
        plotter.add_text(caption, position=(24, 30), font_size=13, color='#43566a')

    def finish():
        plotter.view_xy()
        plotter.enable_parallel_projection()
        plotter.camera.position = (2, 0, 12)
        plotter.camera.focal_point = (2, 0, 0)
        plotter.camera.up = (0, 1, 0)
        plotter.camera.parallel_scale = 2.6

    panel(0, 'A  Free-space Gaussian beam',
          'Complex field, Gouy phase, power quadrature\nand coherent beamlet convergence')
    x = np.linspace(0, 4, 200)
    width = .17 * np.sqrt(1 + ((x - 1) / 1.05) ** 2)
    line([(0, 0), (4.3, 0)], grey, 2)
    for sign in (-1, 1):
        line(np.column_stack((x, sign * width)))
    line([(1, -.8), (1, .8)], grey, 2)
    label((.65, .9), 'waist w0 = 0.1 m')
    line([(3.1, -.9), (3.1, .9)], green, 3)
    label((2.6, 1.05), 'receiver samples', green)
    label((2.9, -.95), 'propagation +z')
    finish()

    panel(1, 'B  Plane reflections',
          'Direct 1/r and image-path spreading;\none- and two-plane specular paths')
    line([(2.8, -1.25), (2.8, 1.25)], grey, 8)
    line([(.3, -.7), (2.8, 0)])
    line([(2.8, 0), (.8, .56)], orange)
    dot((.3, -.7)); dot((.8, .56), orange)
    label((.15, -.98), 'Tx')
    label((.55, .8), 'Rx', orange)
    label((2.95, .75), 'PEC plane')
    label((1.3, -.65), 'a')
    label((1.65, .5), 'b')
    finish()

    panel(2, 'C  Spherical mirrors',
          'Normal/oblique, convex/concave reflection;\nastigmatism, caustic phase, repeated foci')
    t = np.linspace(-.67, .67, 150)
    arc = np.column_stack((1.4 + 2 * np.cos(t), 2 * np.sin(t)))
    line(arc, grey, 7)
    for y in (-.5, .5):
        hit = (1.4 + np.sqrt(4 - y*y), y)
        line([(.1, y), hit])
        line([hit, (1.3, -y * 1.15)], orange)
    dot((2.4, 0), orange)
    label((1.7, -.25), 'paraxial focus', orange)
    label((.3, 1.1), 'concave section shown')
    label((2.9, -1.45), 'R = 2 m')
    finish()

    panel(3, 'D  Cylindrical mirrors',
          'One curved transverse direction;\none-axis caustics and repeated returns')
    line(arc, grey, 7)
    for y in (-.5, .5):
        hit = (1.4 + np.sqrt(4 - y*y), y)
        line([(.1, y), hit])
        line([hit, (1.3, -y * 1.15)], orange)
    label((.1, 1.1), 'section normal to cylinder axis')
    label((.25, -1.15), 'axis out of page: no focusing')
    label((2.95, -1.5), 'R = 2 m')
    finish()

    panel(4, 'E  Dielectric slab',
          'Coherent Fabry-Perot sum, loss, Fresnel/Snell;\ninteraction-depth and branch-pruning sweeps')
    plotter.add_mesh(pv.Plane(center=(2, 0, -.05), direction=(0, 0, 1),
                              i_size=1.2, j_size=2.5), color='#d3eae5', lighting=False)
    for boundary in (1.4, 2.6):
        line([(boundary, -1.25), (boundary, 1.25)], green, 3)
    line([(.1, .5), (1.4, .5)])
    line([(1.4, .5), (2.6, .5), (1.4, 0), (2.6, -.5)], orange)
    line([(2.6, .5), (4.25, .5)], green)
    line([(2.6, -.5), (4.25, -.5)], green)
    label((.1, .8), 'incident')
    label((1.55, .9), 'epsilon_r = 4')
    label((1.45, -1.5), 'thickness 0.013 m')
    label((2.85, .8), 'transmitted', green)
    finish()

    panel(5, 'F  Curved dielectric interface',
          'Gaussian refraction: complex field and width;\nn = 1.25, 1.5, 2.0 against ABCD optics')
    t = np.linspace(np.pi - .65, np.pi + .65, 150)
    line(np.column_stack((3.6 + 2 * np.cos(t), 2 * np.sin(t))), green, 7)
    for y in (-.5, .5):
        hit = (3.6 - np.sqrt(4 - y*y), y)
        line([(.05, y), hit])
        line([hit, (4.2, y * .4)], green)
    label((.1, .9), 'vacuum')
    label((2.25, 1.15), 'dielectric sphere')
    label((2.35, -.95), 'R = 2 m')
    finish()
    target = output / 'geometries.png'
    plotter.screenshot(str(target))
    plotter.close()
    print(target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'figures')
    render(parser.parse_args().output)
