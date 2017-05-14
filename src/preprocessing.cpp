/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file image2mesh
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2014/04/22
 *
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//#include "DGTal/io/writer/GenericWriter.h"
//#include "DGTal/io/readers/GenericWriter.h"

#include <QtGui/qapplication.h>
#include <DGtal/io/writers/PGMWriter.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include "DGtal/io/colormaps/HueShadeColorMap.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

typedef ImageContainerBySTLVector< Z2i::Domain, unsigned char> Image2D;

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info()<<std::endl;
  exit(1);
}


void erosion(Image2D & image) {
  // struct elt: 3x3 block


  Z2i::Point upperBound = image.domain().upperBound();


  for (int i = 1; i < upperBound[0] - 1; i++)
    for (int j = 1; j < upperBound[1] - 1; j++) 
      if (  (int)image( Z2i::Point(i-1, j-1)) > 0 
         && (int)image( Z2i::Point(i-1, j  )) > 0
         && (int)image( Z2i::Point(i-1, j+1)) > 0
         && (int)image( Z2i::Point(i,   j-1)) > 0
         && (int)image( Z2i::Point(i,   j  )) > 0
         && (int)image( Z2i::Point(i,   j+1)) > 0
         && (int)image( Z2i::Point(i+1, j-1)) > 0
         && (int)image( Z2i::Point(i+1, j  )) > 0
         && (int)image( Z2i::Point(i+1, j+1)) > 0
         )
        image.setValue(Z2i::Point(i, j), 2);

  for (int i = 1; i < upperBound[0] - 1; i++)
    for (int j = 1; j < upperBound[1] - 1; j++) 
      if (  (int)image( Z2i::Point(i, j)) == 2) 
        image.setValue(Z2i::Point(i, j), 255);
      else 
        image.setValue(Z2i::Point(i, j), 0);

}

void dilation(Image2D & image) {
  // struct elt: 3x3 block


  Z2i::Point upperBound = image.domain().upperBound();


  for (int i = 1; i < upperBound[0] - 1; i++)
    for (int j = 1; j < upperBound[1] - 1; j++) 
      if (  (int)image( Z2i::Point(i, j)) == 255 ) {
        image.setValue(Z2i::Point(i-1, j-1), 2);
        image.setValue(Z2i::Point(i-1, j  ), 2);
        image.setValue(Z2i::Point(i-1, j+1), 2);
        image.setValue(Z2i::Point(i,   j-1), 2);
        image.setValue(Z2i::Point(i,   j  ), 2);
        image.setValue(Z2i::Point(i,   j+1), 2);
        image.setValue(Z2i::Point(i+1, j-1), 2);
        image.setValue(Z2i::Point(i+1, j  ), 2);
        image.setValue(Z2i::Point(i+1, j+1), 2);
      }

  for (int i = 1; i < upperBound[0] - 1; i++)
    for (int j = 1; j < upperBound[1] - 1; j++) 
      if (  (int)image( Z2i::Point(i, j)) == 2 )
        image.setValue(Z2i::Point(i, j), 255);
      else 
        image.setValue(Z2i::Point(i, j), 0);

}

inline void opening(Image2D & image) {
  erosion(image);
  dilation(image);
}

inline void closing(Image2D & image) {
  dilation(image);
  erosion(image);
}

int main(int argc, char **argv)
{
  if (argc == 1)
    missingParam("path to image");

  //Reading the image
  typedef ImageContainerBySTLVector< Z2i::Domain, unsigned char> Image;
  Image2D image = GenericReader<Image>::import(argv[1]);

  opening(image);
  closing(image);

  //GenericWriter<Image>::export("outputtest", image);
  image >> "clo(ope).pgm";

  image = GenericReader<Image>::import(argv[1]);
  closing(image);
  opening(image);
  image >> "ope(clo).pgm";

  image = GenericReader<Image>::import(argv[1]);
  closing(image);
  image >> "clo.pgm";

  image = GenericReader<Image>::import(argv[1]);
  opening(image);

  image >> "ope.pgm";
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
