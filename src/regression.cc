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
// along with regression.  If not, see <https://www.gnu.org/licenses/>.
//
// Copyright (C) 2020 Greg Hedger
//

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <vector>

using namespace std;

namespace hedger {
  struct DataPoint {
    double x;
    double y;
  };

  void printUsage() {
    printf("regression\n");
    printf("Ordinary Least Squares (OLS) linear regression analysis.\n");
    printf("Calculates Y baseline b and slope m from set of {x,y} points.\n");
    printf("Copyright (C) 2020 Greg Hedger\n");
    printf("\nOptions:\n");
    printf("  -f Specify CSV or other non-digit-separated file\n");
    printf("  -xf Specify file and swap x and y values\n");
    printf("\nUsage:\n");
    printf(" regression [x₁] [y₁] ... [xₙ] [yₙ]\n");
    printf(" regression -f [csv_file]\n");
    printf(" regression -xf [csv_file]\n");
    printf("CSV files can use any non-digit separator.");
  }

  // getSums
  // Get requisite sums (sigmas) for the best fit calculations
  // Entry: data array of {x,y} points
  //        size of array
  //        pointer to destination x variable
  //        pointer to destination y variable
  //        pointer to sum of x²
  //        pointer to sum of xy
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

  // getMean
  // Get x mean (x̄)
  // Entry: data array of {x,y} points
  //        size of array
  double getMean(
      hedger::DataPoint *data,
      size_t size
      )
  {
    double sum = 0.0;
    for( size_t i = 0; i < size; i++ ) {
      sum += data[i].x;
    }
    sum /= size;
    return sum;
  }

  // getBestFit
  //
  // Get the best slope m and baseline b from a set of {x,y} points.
  //
  // Entry: pointer to data
  //        # of datum
  //        pointer to baseline b result
  //        pointer to slope m result
  void getBestFit(hedger::DataPoint *data, size_t size, double *b, double *m)
  {
    double sigmaX = 0.0, sigmaY = 0.0;
    double sigmaXY = 0.0;
    double sigmaXSquared = 0.0;
    // Sum x and y and their squares
    getSums(data, size, &sigmaX, &sigmaY, &sigmaXSquared, &sigmaXY);

    //
    //         N Σ(xy) − Σx Σy
    //    m = -----------------
    //         N Σ(x²) − (Σx)²
    //
    //         Σy − m Σx
    //    b = -----------
    //            N
    //
    //
    //    y = mx + b
    //
    *m = size * (sigmaXY) - (sigmaX * sigmaY);
    *m /= (size * sigmaXSquared) - sigmaX * sigmaX;
    *b = sigmaY - *m * sigmaX;
    *b /= size;
  }

  // getLeastSquares
  //
  // Get ordinary least squares regression on a set of { x, y } data
  // (slight rearrangement of getBestFit equation)
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
    getSums(data, size, &sigmaX, &sigmaY, &sigmaXSquared, &sigmaXY);

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

  // parseFile
  // Parses a CSV or other non-numeric value separated file and gets the data into an array of DataPoints.
  // Entry: filename
  // Exit: pointer to DataPoint array
  DataPoint *parseFile(const char *file, int *data_size) {
    vector<double> x, y;
    int tuples = 0;
    bool xy = false;
    const int MAX_DIGITS = 256;
    char accum[MAX_DIGITS], c;
    int aidx = 0; // accumulator index
    double d = 0.0;
    DataPoint *data = NULL;
    FILE *f = fopen(file, "r");
    if (f) {
      while ((c = fgetc(f)) != EOF) {
        //printf("%c", c);
        if (isdigit(c) || '.' == c || '-' == c) {
          if (aidx > MAX_DIGITS) {
            x.clear();
            y.clear();
            // Error;
            return NULL;
          }
          accum[aidx++] = c; //  read into number accumulator
        } else {
          // Nondigit,nunumeric; terminate, read whole # from accumulator
          accum[aidx] = '\0'; // null terminate
          aidx = 0; // reset accum index
          sscanf(accum, "%lf", &d);
          if (xy) {
            y.push_back(d);
            tuples++;
          } else {
            x.push_back(d);
          }
          xy ^= true; // Toggle x/y
        }
      }
      // Allocate data source and free up lists
      data = new DataPoint[tuples];
      for (auto i = 0; i < tuples; i++) {
        data[i].x = x[i];
        data[i].y = y[i];
      }
      x.clear();
      y.clear();
      // Set data size
      *data_size = tuples;
    }
    return data;
  }

} // namespace hedger

static const int DATA_SIZE = 6;
int main(int argc, const char *argv[])
{
  using namespace hedger;
  DataPoint *data = NULL;
  int data_size = 0;

  if (argc < 5) {
    if (argc > 2) {
      if (!strcmp(argv[1], "-f") || !strcmp(argv[1], "-xf")) {
        data = parseFile(argv[2], &data_size);
        if (!data) {
          printf("Could not read or allocate data, file '%s'\n", argv[2]);
          return -1;
        }
        if (!strcmp(argv[1], "-xf")) {
          for (auto i = 0; i < data_size; i++) {
            double swap = data[i].y;
            data[i].y = data[i].x;
            data[i].x = swap;
          }
        }
      }
    } else {
      printUsage();
      return 1;
    }
  }

  // Pull in arguments from command line, if and only if we are not
  // parsing from a file.
  if (!data) {
    data_size = (argc - 1) / 2;
    data = new DataPoint[ data_size ];
    double v;
    for (int i = 0; i < data_size; i++) {
      sscanf(argv[ 1 + i * 2 ], "%lf", &v);
      data[ i ].x = v;
      sscanf(argv[ 1 + i * 2 + 1 ], "%lf", &v);
      data[ i ].y = v;
    }
  }

  // Warn if odd # of arguments
  if (data_size * 2 < (argc - 1)) {
    printf("\nWARNING: Ignoring last param!");
  }

  // Get the Y baseline ("b") and slope ("m")
  double m = 0.0, b = 0.0;
  getBestFit( data, data_size, &b, &m );
  printf("Best fit (OLS):\n");
  printf("b=%lf\nm=%lf\n", b, m);

  double xbar = getMean( data, data_size );

  // Print y at center point x-bar
  printf("\ny=%lf at x=x̄=%lf\n", m*xbar + b, xbar);

  // Free resources and exit
  delete data;
  return 0;
}
