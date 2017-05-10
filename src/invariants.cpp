///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/readers/PGMReader.h>
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
///////////////////////////////////////////////////////////////////////////////

#define SQ(x) pow(x, 2)

using namespace std;
using namespace DGtal;
using namespace DGtal::Z2i; //We'll only consider Z² digital space on
				//32bit integers
typedef array<double, 2> double2;

//Image type (image of unsigned int)
typedef ImageContainerBySTLVector< Domain, unsigned int > Image;

double ccw(const Point &a, const Point &b, const Point &c)
{
	Point u = b - a, v = c - a;
	return u[0]*v[1] - u[1]*v[0];
}

double M(const DigitalSet &S, int p, int q, Point c)
{
	double moment = 0;
	for(auto it = S.begin(), itEnd = S.end(); it != itEnd; ++it)
		moment += pow((*it)[0]-c[0], p) * pow((*it)[1]-c[1], q);
	return moment;
}

bool ACTIVATE_SCALE = 0;
double M0(const DigitalSet &S, int p, int q) { return M(S, p, q, Point(0, 0));   }
Point massCenter(const DigitalSet &S)
{
	return Point(M0(S,1,0) / S.size(), M0(S,0,1) / S.size());
}
double M (const DigitalSet &S, int p, int q) { return M(S, p, q, massCenter(S)) / pow(S.size(), ACTIVATE_SCALE*(1+p/2+q/2)); }





vector<double> invariants(const DigitalSet &S)
{
	double m11 = M(S, 1, 1);
	double m20 = M(S, 2, 0),    m02 = M(S, 0, 2);
	double m30 = M(S, 3, 0),    m03 = M(S, 0, 3);
	double m12 = M(S, 1, 2),    m21 = M(S, 2, 1);
	
	double m2103 = m21+m03,     s2103 = SQ(m2103);
	double m3012 = m30+m12,     s3012 = SQ(m3012);
	double m_3 = m30-3*m12,     b1 = m3012*(s3012 - 3*s2103);
	double m3_ = 3*m21-m03,     b2 = m2103*(3*s3012 - s2103);

	double f1 = m20 + m02;
	double f2 = SQ(m20-m02) + 4*SQ(m11);
	double f3 = SQ(m_3) + SQ(m3_);
	double f4 = s3012 + s2103;
	double f6 = (m20-m02)*(s3012-s2103) + 4*m11*m3012*m2103;
	double f5 = m_3*b1 + m3_*b2;
	double f7 = m3_*b1 - m_3*b2;

	return vector<double>{f1, f2, f3, f4, f5, f6, f7};
}


// in trouble w/ noise
// invariant by resize / rotation
// helps estimating the thickness of the object
vector<double> rappPerimAire (const Image &image) {
	double perim (0), aire (0);

  Z2i::Point upperBound = image.domain().upperBound();

	double im_size (upperBound[0]*upperBound[1]);
	
	for (int i = 1; i < upperBound[0] -1; i++) {
		for (int j = 1; j < upperBound[1] -1; j++) {
			if ((int)image( Z2i::Point(i-1, j-1)) != 0) {
				aire++;

	      if (  (int)image( Z2i::Point(i-1, j-1)) == 0 
	         || (int)image( Z2i::Point(i-1, j  )) == 0
	         || (int)image( Z2i::Point(i-1, j+1)) == 0
	         || (int)image( Z2i::Point(i,   j-1)) == 0
	         || (int)image( Z2i::Point(i,   j+1)) == 0
	         || (int)image( Z2i::Point(i+1, j-1)) == 0
	         || (int)image( Z2i::Point(i+1, j  )) == 0
	         || (int)image( Z2i::Point(i+1, j+1)) == 0
	         )
	      	perim++;
			}
		}  
	}
	return vector<double>{(perim/aire), (aire/im_size)};
}

// bit of cheating: representing the % of area covered by the object in the image

void invariants_image(const Image &image, vector<double> &inv) {
	vector<double> temp (rappPerimAire(image));
	inv.push_back(temp[0]);
	inv.push_back(temp[1]);
}

int main(int argc, char *argv[])
{

	string filename = (argc < 2) ? "../database/bat-5.pgm" : argv[1];


	//We read the PGM file
	Image image = PGMReader<Image>::importPGM(filename);
	//trace.info() << "Image read :"<< image << endl;

	//We convert pixels in ]0,255] into a digital set
	DigitalSet set2d( image.domain() );
	SetFromImage<DigitalSet>::append<Image>(set2d, image, 0, 255);

	vector<double> inv = invariants(set2d);
    invariants_image(image, inv);
	cout << inv[0];
	for(int i=1; i<inv.size(); i++)
	{
		cout << "," << inv[i];
	}

	return 0;
}
