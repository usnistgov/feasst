/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include <limits.h>
#include "base_random.h"
#ifdef GSL_
  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <gsl/gsl_errno.h>
  #include <gsl/gsl_spline.h>
#endif  // GSL_

using namespace feasst;

TEST(Base, base) {
  Base b;
}

#ifdef JSON_
TEST(json, json) {
  // create an empty structure (null)
 nlohmann::json j;

 // add a number that is stored as double (note the implicit conversion of j to an object)
 j["pi"] = 3.141;

 // add a Boolean that is stored as bool
 j["happy"] = true;

 // add a string that is stored as std::string
 j["name"] = "Niels";

 // add another null object by passing nullptr
 j["nothing"] = nullptr;

 // add an object inside the object
 j["answer"]["everything"] = 42;
 j["answer"]["nothing"] = "everything";

 // add an array that is stored as std::vector (using an initializer list)
 j["list"] = { 1, 0, 2 };

 // add another object (using an initializer list of pairs)
 j["object"] = { {"currency", "USD"}, {"value", 42.99} };

  //// print json with given number of spaces for indent
  //std::cout << j.dump(2) << std::endl;
  //cout << std::setw(4) << j << endl;

  // write prettified JSON to another file
  std::ofstream o("tmp/pretty.json");
  o << std::setw(4) << j << std::endl;

  // read a JSON file
  std::ifstream i("tmp/pretty.json");
  nlohmann::json j2;
  i >> j2;

  // write prettified JSON to another file
  std::ofstream o2("tmp/pretty2.json");
  o2 << std::setw(2) << j2 << std::endl;

  EXPECT_EQ(j["pi"], 3.141);
  EXPECT_TRUE(j["name"] == "Niels");
  EXPECT_EQ(int(j["list"].size()), 3);
  EXPECT_EQ(int(j["object"].size()), 2);
  EXPECT_TRUE(j["object"]["value"] == 42.99);
  EXPECT_TRUE(j["object"]["currency"] == "USD");

  //for (int i = 0; i < j["answer"].size(); ++i) {
  //  cout << j["answer"][i] << endl;
  //}
  nolohmann::json jset = j2["answer"];
  std::ofstream o3("tmp/pretty3.json");
  o3 << std::setw(2) << jset << std::endl;
  EXPECT_TRUE(jset["nothing"] == "everything");

}
#endif  // JSON_

#ifdef HDF5_
TEST(hdf5, hdf5) {
  string fileName("tmp/test.hdf5");

  /*
   * Create a file.
   */
  const H5std_string FILE_NAME( "tmp/Select.h5" );
  H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

  /*
  * Create property list for a dataset and set up fill values.
  */
  int fillvalue = 0;   /* Fill value for the dataset */
  DSetCreatPropList plist;
  plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

  /*
   * Create dataspace for the dataset in the file.
   */
  const int   FSPACE_RANK = 2;  // Dataset rank as it is stored in the file
  const int   FSPACE_DIM1 = 8;  // Dimension sizes of the dataset as it is
  const int   FSPACE_DIM2 = 12;  //   stored in the file
  hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2}; // dim sizes of ds (on disk)
  DataSpace fspace( FSPACE_RANK, fdim );

  /*
   * Create dataset and write it into the file.
   */
  const H5std_string DATASET_NAME( "Matrix in file" );
  DataSet* dataset = new DataSet(file->createDataSet(
    DATASET_NAME, PredType::NATIVE_INT, fspace, plist));

  /*
   * Select hyperslab for the dataset in the file, using 3x2 blocks,
   * (4,3) stride and (2,4) count starting at the position (0,1).
   */
  hsize_t start[2]; // Start of hyperslab
  hsize_t stride[2]; // Stride of hyperslab
  hsize_t count[2];  // Block count
  hsize_t block[2];  // Block sizes
  start[0]  = 0; start[1]  = 1;
  stride[0] = 4; stride[1] = 3;
  count[0]  = 2; count[1]  = 4;
  block[0]  = 3; block[1]  = 2;
  fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

  /*
   * Create dataspace for the first dataset.
   */
  const int   MSPACE1_RANK = 1;  // Rank of the first dataset in memory
  const int   MSPACE1_DIM = 50;   // Dataset size in memory
  hsize_t dim1[] = {MSPACE1_DIM};  /* Dimension size of the first dataset
                                     (in memory) */
  DataSpace mspace1( MSPACE1_RANK, dim1 );

  /*
   * Select hyperslab.
   * We will use 48 elements of the vector buffer starting at the
   * second element.  Selected elements are 1 2 3 . . . 48
   */
  start[0]  = 1;
  stride[0] = 1;
  count[0]  = 48;
  block[0]  = 1;
  mspace1.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

  /*
   * Write selection from the vector buffer to the dataset in the file.
   *
   * File dataset should look like this:
   *                    0  1  2  0  3  4  0  5  6  0  7  8
   *                    0  9 10  0 11 12  0 13 14  0 15 16
   *                    0 17 18  0 19 20  0 21 22  0 23 24
   *                    0  0  0  0  0  0  0  0  0  0  0  0
   *                    0 25 26  0 27 28  0 29 30  0 31 32
   *                    0 33 34  0 35 36  0 37 38  0 39 40
   *                    0 41 42  0 43 44  0 45 46  0 47 48
   *                    0  0  0  0  0  0  0  0  0  0  0  0
   */
  int    vector[MSPACE1_DIM];  // vector buffer for dset

  /*
   * Buffer initialization.
   */
  vector[0] = vector[MSPACE1_DIM - 1] = -1;
  for (int i = 1; i < MSPACE1_DIM - 1; i++)
      vector[i] = i;

  dataset->write( vector, PredType::NATIVE_INT, mspace1, fspace );

  /*
   * Reset the selection for the file dataspace fid.
   */
  fspace.selectNone();

  /*
   * Create dataspace for the second dataset.
   */
  const int   MSPACE2_RANK = 1;  // Rank of the second dataset in memory
  const int   MSPACE2_DIM = 4;  // Dataset size in memory
  hsize_t dim2[] = {MSPACE2_DIM};  /* Dimension size of the second dataset
                                     (in memory */
  DataSpace mspace2( MSPACE2_RANK, dim2 );

  /*
   * Select sequence of NPOINTS points in the file dataspace.
   */
  const int   NPOINTS = 4;  // Number of points that will be selected
  hsize_t coord[NPOINTS][FSPACE_RANK]; /* Array to store selected points
                                          from the file dataspace */
  coord[0][0] = 0; coord[0][1] = 0;
  coord[1][0] = 3; coord[1][1] = 3;
  coord[2][0] = 3; coord[2][1] = 5;
  coord[3][0] = 5; coord[3][1] = 6;

  fspace.selectElements( H5S_SELECT_SET, NPOINTS, (const hsize_t *)coord);

  /*
   * Write new selection of points to the dataset.
   */
  int    values[] = {53, 59, 61, 67};  /* New values to be written */
  dataset->write( values, PredType::NATIVE_INT, mspace2, fspace );

  /*
   * File dataset should look like this:
   *                   53  1  2  0  3  4  0  5  6  0  7  8
   *                    0  9 10  0 11 12  0 13 14  0 15 16
   *                    0 17 18  0 19 20  0 21 22  0 23 24
   *                    0  0  0 59  0 61  0  0  0  0  0  0
   *                    0 25 26  0 27 28  0 29 30  0 31 32
   *                    0 33 34  0 35 36 67 37 38  0 39 40
   *                    0 41 42  0 43 44  0 45 46  0 47 48
   *                    0  0  0  0  0  0  0  0  0  0  0  0
   *
   */

  /*
   * Close the dataset and the file.
   */
  delete dataset;
  delete file;

}
#endif  //HDF5_

#ifdef HDF5_
TEST(hdf5, compress) {
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                 *
 * Copyright by the Board of Trustees of the University of Illinois.       *
 * All rights reserved.                   *
 *                                                                       *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have       *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 *  This example illustrates how to create a compressed dataset.
 *  It is used in the HDF5 Tutorial.
 */

const H5std_string  FILE_NAME("tmp/h5tutr_cmprss.h5");
const H5std_string  DATASET_NAME("Compressed_Data");
const int  DIM0 = 100;
const int  DIM1 = 20;

    hsize_t dims[2] = { DIM0, DIM1 };  // dataset dimensions
    hsize_t chunk_dims[2] = { 20, 20 };  // chunk dimensions
    int     i,j, buf[DIM0][DIM1];

    // Try block to detect exceptions raised by any of the calls inside it
    try
    {
  // Turn off the auto-printing when failure occurs so that we can
  // handle the errors appropriately
  Exception::dontPrint();

  // Create a new file using the default property lists.
  H5File file(FILE_NAME, H5F_ACC_TRUNC);

  // Create the data space for the dataset.
  DataSpace *dataspace = new DataSpace(2, dims);

  // Modify dataset creation property to enable chunking
  DSetCreatPropList  *plist = new  DSetCreatPropList;
  plist->setChunk(2, chunk_dims);

  // Set ZLIB (DEFLATE) Compression using level 6.
  // To use SZIP compression comment out this line.
  plist->setDeflate(6);

  // Uncomment these lines to set SZIP Compression
  // unsigned szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  // unsigned szip_pixels_per_block = 16;
  // plist->setSzip(szip_options_mask, szip_pixels_per_block);

  // Create the dataset.
  DataSet *dataset = new DataSet(file.createDataSet( DATASET_NAME,
                          PredType::STD_I32BE, *dataspace, *plist) );

  for (i = 0; i< DIM0; i++)
    for (j=0; j<DIM1; j++)
        buf[i][j] = i+j;

  // Write data to dataset.
  dataset->write(buf, PredType::NATIVE_INT);

  // Close objects and file.  Either approach will close the HDF5 item.
  delete dataspace;
  delete dataset;
  delete plist;
  file.close();

  // -----------------------------------------------
  // Re-open the file and dataset, retrieve filter
  // information for dataset and read the data back.
  // -----------------------------------------------

  int        rbuf[DIM0][DIM1];
  int        numfilt;
  size_t     nelmts={1}, namelen={1};
  unsigned  flags, filter_info, cd_values[1], idx;
  char       name[1];
  H5Z_filter_t filter_type;

  // Open the file and the dataset in the file.
  file.openFile(FILE_NAME, H5F_ACC_RDONLY);
  dataset = new DataSet(file.openDataSet( DATASET_NAME));

  // Get the create property list of the dataset.
  plist = new DSetCreatPropList(dataset->getCreatePlist ());

  // Get the number of filters associated with the dataset.
  numfilt = plist->getNfilters();
  cout << "Number of filters associated with dataset: " << numfilt << endl;

  for (idx=0; int(idx) < numfilt; idx++) {
      nelmts = 0;

      filter_type = plist->getFilter(idx, flags, nelmts, cd_values, namelen, name , filter_info);

      cout << "Filter Type: ";

      switch (filter_type) {
        case H5Z_FILTER_DEFLATE:
             cout << "H5Z_FILTER_DEFLATE" << endl;
             break;
        case H5Z_FILTER_SZIP:
             cout << "H5Z_FILTER_SZIP" << endl;
             break;
        default:
             cout << "Other filter type included." << endl;
        }
  }

  // Read data.
  dataset->read(rbuf, PredType::NATIVE_INT);

  delete plist;
  delete dataset;
  file.close();  // can be skipped

    }  // end of try block

    // catch failure caused by the H5File operations
    catch(FileIException error)
    {
  error.printError();
  //return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
  error.printError();
  //return -1;
    }

    // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error)
    {
  error.printError();
  //return -1;
    }

    //return 0;  // successfully terminated
}
#endif  //HDF5_

#ifdef GSL_
TEST(Base, GSL) {

  int i;
  double x[10], y[10];
  //double xi, yi, x[10], y[10];

  //printf ("#m=0,S=2\n");

  for (i = 0; i < 10; i++)
    {
      x[i] = i + 0.5 * sin (i);
      y[i] = i + cos (i * i);
      //printf ("%g %g\n", x[i], y[i]);
    }

  //printf ("#m=1,S=0\n");

  {
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, 10);

    gsl_spline_init (spline, x, y, 10);

//    for (xi = x[0]; xi < x[9]; xi += 0.01) {
//      yi = gsl_spline_eval (spline, xi, acc);
//      //printf ("%g %g\n", xi, yi);
//    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
}
#endif  //GSL_

