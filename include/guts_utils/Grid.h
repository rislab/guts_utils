#pragma once

/*
Copyright 2024 Wennie Tabib

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

// Implementation of
// [1] Amanatides, John, and Andrew Woo. "A fast voxel traversal
// algorithm for ray tracing." Eurographics. Vol. 87. No. 3. 1987.

namespace guts_utils
{

  typedef struct Cell
  {
    unsigned int row;
    unsigned int col;
    Cell(unsigned int r, unsigned int c)
    {
      row = r;
      col = c;
    }

  } cell_t;

  typedef struct Grid
  {

    typedef std::shared_ptr<Grid> Ptr;
    typedef std::shared_ptr<const Grid> ConstPtr;
    typedef Eigen::Vector2f Point;

    unsigned int width;
    unsigned int height;
    float resolution;
    Point origin;
    std::vector<float> data;

    Grid() {}

    Grid(Point orig, unsigned int w, unsigned int h, float res)
    {
      width = w;
      height = h;
      resolution = res;
      origin = orig;
      data.resize( w * h, 0.0);
    }

    void set(const Cell& c, float value)
    {
      unsigned int index = cellToIndex(c);

      if (index < width*height)
	data[index] = value;
      else
	std::cerr << "Error! grid is smaller (" << width << "x" << height <<
	  " than requested index: " << index<< std::endl;
    }

    Cell pointToCell(const Point& p) const
    {
      return Cell((unsigned int)((1.0f/resolution)*(p.y() - origin.y())),
		  (unsigned int)((1.0f/resolution)*(p.x() - origin.x())));
    }

    unsigned int pointToIndex(const Point& p) const
    {
      return cellToIndex( pointToCell(p) );
    }

    Point cellToPoint(const Cell& c) const
    {
      return Point( (c.col * resolution) + origin.x(),
		    (c.row * resolution) + origin.y() );
    }

    unsigned int cellToIndex(const Cell& c) const
    {
      return c.row * width + c.col;
    }

    Cell indexToCell(const unsigned int index) const
    {
      return Cell(index / width, index & width);
    }

    bool inGrid(Point p) const
    {
      Cell c = pointToCell(p);
      return inGrid(c);
    }

    bool inGrid(Cell c) const
    {
      if (c.row < height and c.col < width)
	return true;
      return false;
    }

    bool inGrid(unsigned int index) const
    {
      Cell c = indexToCell(index);
      return inGrid(c);
    }


    void update_row (Cell& c, Point& tm, Point& td, int step) const
    {
      c.row += step;
      tm.y() += td.y();
    }

    void update_col (Cell& c, Point& tm, Point& td, int step) const
    {
      c.col += step;
      tm.x() += td.x();
    };

    std::pair<bool, std::vector<Cell>> raytrace(const Point& start, const Point& end) const
    {
      std::vector<Cell> traversed_cells;

      // Check if ray origin is in the grid
      if (!inGrid(start))
      {
	std::cerr << "Ray origin does not lie within the grid" << std::endl;
	return std::make_pair(false, traversed_cells);
      }

      if (!inGrid(end))
      {
	std::cerr << "Ray endpoint does not lie within the grid" << std::endl;
	return std::make_pair(false, traversed_cells);
      }

      // Determine the cell in which the ray origin is found
      Cell corigin = pointToCell(start);
      Cell cend = pointToCell(end);

      // If start and end cell are the same, return start
      if ( (corigin.row == cend.row) && (corigin.col == cend.col) )
      {
	traversed_cells.push_back(corigin);
	return std::make_pair(true, traversed_cells);
      }

      // traversed_cells should contain more than one cell, so we
      // can reserve some number of cells
      unsigned int num_traversed_cells = std::abs( ((int)corigin.row) - ((int)cend.row) )
	+ std::abs( ((int)corigin.col) - ((int)cend.col) ) + 1;
      traversed_cells.reserve(num_traversed_cells);

      if (num_traversed_cells > width * height)
      {
	std::cerr << "Total number of cells larger than contained in grid." << std::endl;
	return std::make_pair(false, traversed_cells);
      }

      Point direction = (end - start).normalized();
      Point tMax;
      Point tDelta;

      Point c = cellToPoint(corigin);

      // Check if we can restrict our search to be along row or along column
      bool search_along_col = !(corigin.col == cend.col);
      bool search_along_row = !(corigin.row == cend.row);

      // Initialize step_row (y) and step_col (x) to be either 1 or -1
      // to indicate whether to cross cell boundaries as determined by
      // the sign of the x- and y-components of the direction vector.

      // From [1]:
      // step_col == stepX
      int step_col = -1;
      if (direction.x() > 0.0f)
      {
	step_col = 1;
	c.x() += resolution;
      }

      // From [1]:
      // step_row == stepY
      int step_row = -1;
      if (direction.y() > 0.0f)
      {
	step_row = 1;
	c.y() += resolution;
      }

      if (search_along_col)
      {
	tDelta.x() = ( resolution * ((float) step_col) ) / direction.x();
	tMax.x() = (c.x() - start.x()) / direction.x();
      }

      if (search_along_row)
      {
	tDelta.y() = ( resolution * ((float) step_row) ) / direction.y();
	tMax.y() = (c.y() - start.y()) / direction.y();
      }

      Cell currc = corigin;
      while (true)
      {
	traversed_cells.push_back(currc);

	if ( (currc.row == cend.row) && (currc.col == cend.col) )
	  return std::make_pair(true, traversed_cells);

	if ( search_along_row == true && search_along_col == false)
	  update_row(currc, tMax, tDelta, step_row);

	else if ( search_along_row == false && search_along_col == true)
	  update_col(currc, tMax, tDelta, step_col);

	else if (search_along_row == true && search_along_col == true)
	{
	  if (tMax.x() < tMax.y())
	    update_col(currc, tMax, tDelta, step_col);
	  else
	    update_row(currc, tMax, tDelta, step_row);
	}

	else // this case should never happen
	{
	  std::cerr << "Raytracing hit an impossible case...returning." << std::endl;
	  return std::make_pair(false, traversed_cells);
	}

      }
      return std::make_pair(true, traversed_cells);
    }

  } grid_t;

}
