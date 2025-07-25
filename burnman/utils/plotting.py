# MIT License

# Copyright (c) 2017 Michal Haták
# Copyright (c) 2025 Bob Myhill

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import heapq


def _get_seg_dist_sq(px, py, a, b):
    """
    Compute the squared distance from point (px, py) to the line segment defined by points a and b.

    :param px: X-coordinate of the point.
    :param py: Y-coordinate of the point.
    :param a: First point of the segment as a tuple (x, y).
    :param b: Second point of the segment as a tuple (x, y).
    :return: Squared distance from point to segment.
    :rtype: float
    """
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay

    if dx != 0 or dy != 0:
        t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
        if t > 1:
            ax, ay = bx, by
        elif t > 0:
            ax += dx * t
            ay += dy * t

    return (px - ax) ** 2 + (py - ay) ** 2


def _point_to_polygon_distance(x, y, polygon):
    """
    Compute the signed distance from a point to the nearest polygon edge.

    :param x: X-coordinate of the point.
    :param y: Y-coordinate of the point.
    :param polygon: List of rings (each ring is a numpy array of shape (N, 2)).
    :return: Signed distance; positive if inside, negative if outside.
    :rtype: float
    """
    inside = False
    min_dist_sq = np.inf

    for ring in polygon:
        b = ring[-1]
        for a in ring:
            if ((a[1] > y) != (b[1] > y)) and (
                x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0]
            ):
                inside = not inside
            min_dist_sq = min(min_dist_sq, _get_seg_dist_sq(x, y, a, b))
            b = a

    distance = np.sqrt(min_dist_sq)
    return distance if inside else -distance


class Cell:
    """
    Represents a square cell used during the polygon center search.

    :param x: X-coordinate of the cell center.
    :param y: Y-coordinate of the cell center.
    :param h: Half-size of the cell.
    :param polygon: The input polygon as a list of rings.
    """

    def __init__(self, x, y, h, polygon):
        self.x = x
        self.y = y
        self.h = h
        self.d = _point_to_polygon_distance(x, y, polygon)
        self.max = self.d + h * np.sqrt(2)

    def __lt__(self, other):
        return self.max > other.max


def _get_centroid_cell(polygon):
    """
    Estimate the polygon's centroid as an initial guess.

    :param polygon: List of polygon rings.
    :return: A Cell object located at the estimated centroid.
    :rtype: Cell
    """
    points = polygon[0]
    area = 0.0
    cx = 0.0
    cy = 0.0
    b = points[-1]

    for a in points:
        f = a[0] * b[1] - b[0] * a[1]
        cx += (a[0] + b[0]) * f
        cy += (a[1] + b[1]) * f
        area += f * 3
        b = a

    if area == 0:
        return Cell(points[0][0], points[0][1], 0.0, polygon)

    return Cell(cx / area, cy / area, 0.0, polygon)


def visual_center_of_polygon(polygon_rings, precision=1.0, with_distance=False):
    """
    Compute the pole of inaccessibility (visual center) of a polygon with the specified precision.

    :param polygon_rings: A polygon represented as a list of rings. Each ring is a numpy array of shape (N, 2).
    :param precision: Desired precision. Stops when improvement is less than this value.
    :param with_distance: If True, also return the distance to the closest edge.
    :return: The [x, y] coordinates of the center, and optionally the distance.
    :rtype: list or tuple
    """
    coords = polygon_rings[0]
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("Expected polygon ring to be an Nx2 array")

    min_x, min_y = np.min(coords, axis=0)
    max_x, max_y = np.max(coords, axis=0)

    width = max_x - min_x
    height = max_y - min_y
    cell_size = min(width, height)
    max_dim = max(width, height)

    h = cell_size / 2.0

    # If the cell is much longer than it is wide (or vice-versa),
    # just return the mean of x and y.
    if cell_size < max_dim / 100:
        mean_x = (max_x + min_x) / 2.0
        mean_y = (max_y + min_y) / 2.0
        return ([mean_x, mean_y], 0.0) if with_distance else [mean_x, mean_y]

    # Initialize priority queue
    queue = []
    heapq.heapify(queue)

    x_coords = np.arange(min_x, max_x, cell_size)
    y_coords = np.arange(min_y, max_y, cell_size)

    for x in x_coords:
        for y in y_coords:
            heapq.heappush(queue, Cell(x + h, y + h, h, polygon_rings))

    best_cell = _get_centroid_cell(polygon_rings)

    bbox_cell = Cell(min_x + width / 2, min_y + height / 2, 0.0, polygon_rings)
    if bbox_cell.d > best_cell.d:
        best_cell = bbox_cell

    while queue:
        cell = heapq.heappop(queue)

        if cell.d > best_cell.d:
            best_cell = cell

        if cell.max - best_cell.d <= precision:
            continue

        h = cell.h / 2
        for dx in [-h, h]:
            for dy in [-h, h]:
                heapq.heappush(queue, Cell(cell.x + dx, cell.y + dy, h, polygon_rings))

    result = [best_cell.x, best_cell.y]
    return (result, best_cell.d) if with_distance else result
