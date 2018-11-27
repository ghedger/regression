// regression.cc
//
// This file is part of regression.
//
// Regression is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Sortbench is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with sortbench.  If not, see <https://www.gnu.org/licenses/>.
//

#include <stdio.h>

namespace hedger {
struct DataPoint {
  double x;
  double y;
};
} // namespace hedger

void getSums(
  hedger::DataPoint *data,
  size_t size,
  double *x,
  double *y,
  double *xSquared,
  double *xy
)
{
  *x = 0.0;
  *y = 0.0;
  *xSquared = 0.0;
  *xy = 0.0;
  for( size_t i = 0; i < size; i++ ) {
    *x += data[i].x;
    *y += data[i].y;
    *xSquared += data[i].x * data[i].x;
    *xy += data[i].x * data[i].y;
  }
}

// getLeastSquares
//
// Get ordinary least squares regression on a set of { x, y } data
//
// Entry: pointer to data
//        # of datum
//        pointer to a result
//        pointer to b result
void getLeastSquares(hedger::DataPoint *data, size_t size, double *a, double *b)
{
  double sigmaX = 0.0, sigmaY = 0.0;
  double sigmaXY = 0.0;
  double sigmaXSquared = 0.0;
  // Sum x and y and their squares
  // (TODO: We don't strictly need y^2; optimize out later)
  getSums( data, size, &sigmaX, &sigmaY, &sigmaXSquared, &sigmaXY );

  // Now calculate a and b per standard linear regression equation
  *a = (sigmaY * sigmaXSquared) - (sigmaX * sigmaXY);
  *a /= (size * sigmaXSquared) - (sigmaX * sigmaX);
  *b = (size * sigmaXY) - (sigmaX * sigmaY);
  *b /= (size * sigmaXSquared) - (sigmaX * sigmaX);
}

static const int DATA_SIZE = 6;
int main()
{
  hedger::DataPoint data[DATA_SIZE] =
  { { 43.0, 99.0 },
    { 21.0, 65.0 },
    { 25.0, 79.0 },
    { 42.0, 75.0 },
    { 57.0, 87.0 },
    { 59.0, 81.0 } };

  double a = 0.0, b = 0.0;
  getLeastSquares( data, DATA_SIZE, &a, &b );
  printf("%lf, %lf\n", a, b);
  return 0;
}
