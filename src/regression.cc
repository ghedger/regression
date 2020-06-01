// regression.cc
//
// This file is part of regression.
//
// Regression is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Regression is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with sortbench.  If not, see <https://www.gnu.org/licenses/>.
//
// Copyright (C) 2018 Greg Hedger
//

#include <stdio.h>

namespace hedger {
  struct DataPoint {
    double x;
    double y;
  };

  void printUsage() {
    printf("regression\n");
    printf("Copyright (C) 2020 Greg Hedger\n");
    printf("Usage:\n");
    printf(" regrssion [x1] [y1] [x2] [y2] ...\n");
  }

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
    getSums( data, size, &sigmaX, &sigmaY, &sigmaXSquared, &sigmaXY );

    // Now calculate a and b per standard linear regression equation
    //
    //      (sigma y)(sigma x^2) - (sigma x)(sigma xy)
    // a =  ------------------------------------------
    //          n(sigma x^2) - (sigma x)^2
    //
    //      n(sigma xy) - (sigma x)(sigma y)
    // b =  --------------------------------
    //           n(sigma x^2) - (sigma x)^2
    //
    *a = (sigmaY * sigmaXSquared) - (sigmaX * sigmaXY);
    *a /= (size * sigmaXSquared) - (sigmaX * sigmaX);
    *b = (size * sigmaXY) - (sigmaX * sigmaY);
    *b /= (size * sigmaXSquared) - (sigmaX * sigmaX);
  }
} // namespace hedger

static const int DATA_SIZE = 6;
int main(int argc, const char *argv[])
{
  using namespace hedger;

  if (argc < 2) {
    printUsage();
    return 1;
  }

  // Pull in arguments from command line
  int data_size = (argc - 1) / 2;
  DataPoint *data = new DataPoint[ data_size ];
  double v;
  for (int i = 0; i < data_size; i++) {
    sscanf(argv[ 1 + i * 2 ], "%lf", &v);
    data[ i ].x = v;
    sscanf(argv[ 1 + i * 2 + 1 ], "%lf", &v);
    data[ i ].y = v;
  }

  // Get the Y starting offset ("x") and rise-over-run ("y")
  double a = 0.0, b = 0.0;
  getLeastSquares( data, data_size, &a, &b );
  printf("%lf, %lf\n", a, b);

  // Free resources and exit
  delete data;
  return 0;
}
