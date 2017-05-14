///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/readers/PGMReader.h>
#include <string>

//#include "DGtal/images/imagesSetsUtils/SetFromImage2D.h"
//#include "DGtal/images/Image2DContainerBySTLVector.h"
///////////////////////////////////////////////////////////////////////////////

#define SQ(x) pow(x, 2)

using namespace std;
using namespace DGtal;
using namespace DGtal::Z2i; //We'll only consider Z² digital space on
				//32bit integers
typedef array<double, 2> double2;

//Image2D type (image of unsigned int)
typedef ImageContainerBySTLVector< Domain, unsigned int > Image2D;



//                                                     ____________________________________
//                                                    |                                    |
//                                                    |        preprocessing: noise        |
//                                                    |____________________________________|

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


//                                                     ____________________________________
//                                                    |                                    |
//                                                    |        invariants digitalSet       |
//                                                    |____________________________________|


double ccw(const Z2i::Point &a, const Z2i::Point &b, const Z2i::Point &c)
{
	Z2i::Point u = b - a, v = c - a;
	return u[0]*v[1] - u[1]*v[0];
}

double M(const DigitalSet &S, int p, int q, Z2i::Point c)
{
	double moment = 0;
	for(auto it = S.begin(), itEnd = S.end(); it != itEnd; ++it)
		moment += pow((*it)[0]-c[0], p) * pow((*it)[1]-c[1], q);
	return moment;
}

bool ACTIVATE_SCALE = 0;
double M0(const DigitalSet &S, int p, int q) { return M(S, p, q, Z2i::Point(0, 0));   }
Z2i::Point massCenter(const DigitalSet &S)
{
	return Z2i::Point(M0(S,1,0) / S.size(), M0(S,0,1) / S.size());
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



//                                                     ____________________________________
//                                                    |                                    |
//                                                    |        invariants Image2D            |
//                                                    |____________________________________|

// in trouble w/ noise
// invariant by resize / rotation
// helps estimating the thickness of the object
vector<double> rappPerimAire (const Image2D &image) {
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
	return vector<double>{(perim/aire), (aire/im_size)}; // 2 arg:  bit of cheating: representing the % of area covered by the object in the image

}




//                                   skeleton

// we work in 1-adjacence here
// string :  -------
//           |0|1|2|
//           |7| |3|
//           |6|5|4|
//           -------

// n = number of adjacent cells in the shape
// if n = 1: we always keep our cell
// if n = 2: we erase it if adjacents (8 cases out of 28), accepted
const string accepted2[8] = {"11000000", "01100000", "00110000", "00011000", "00001100", "00000110", "00000011", "10000001"};
// if n = 3: we erase it if they are all adjacents (8 cases too, out of 56)
const string accepted3[8] = {"11100000", "01110000", "00111000", "00011100", "00001110", "00000111", "10000011", "11000001"};
// if n = 4: we erase it if we have a 1x4 (8 cases) or 2x2 situation (12 cases), so 20 cases out of 70)
const string taccepted4[20] = {"11110000", "01111000", "00111100", "00011110", "00001111", "10000111", "11000011", "11100001" // 1x4 cases
					    "11011000", "11001100", "11000110", "01101100", "01100110", "01100011", "00110110", "00110011", "10110001", "00011011", "10011001", "10001101"}; // 2x2 cases
const string accepted4[8] = {"11110000", "01111000", "00111100", "00011110", "00001111", "10000111", "11000011", "11100001"};
const string accepted5[8] = {"11111000", "01111100", "00111110", "00011111", "10001111", "11000111", "11100011", "11110001"};
const string accepted6[8] = {"11111100", "01111110", "00111111", "10011111", "11001111", "11100111", "11110011", "11111001"};

// if n > 4: we intervet bits to check: accepted = rejected

bool isInside(string const accepted[], int size, string const& test) {
	bool b = false;
	for (int i = 0; i < size; i++) {
		if (accepted[i] == test) {
			b = true;
			break;
		}
	}
	return b;
}

void getSkeleton(Image2D &image, vector<Z2i::Point> &skeleton) { // last function on image, we can modify it
	// preprocessing: get all the point of the shape to decrease time complexity:
	for (int i = 1; i < image.domain().upperBound()[0] -1; i++) // there is never a point on the border of the image
		for (int j = 1; j < image.domain().upperBound()[1] -1; j++)
			if ((int)image(Z2i::Point(i, j)) > 0)
				skeleton.push_back(Z2i::Point(i, j));

	bool changeMade = true;
	while (changeMade) { // while we just modified the image
		changeMade = false;
		for (int ind = 0; ind < skeleton.size(); ind++) {
			int nbrNeighbours (0);
			string neighbours ("00000000");
			int i (skeleton[ind][0]), j (skeleton[ind][1]);
      if (  (int)image( Z2i::Point(i-1, j-1)) > 0) {
      	nbrNeighbours ++;
      	neighbours[0] = '1';
      }
      if (  (int)image( Z2i::Point(i-1, j  )) > 0) {
      	nbrNeighbours ++;
      	neighbours[1] = '1';
      }
      if (  (int)image( Z2i::Point(i-1, j+1)) > 0) {
      	nbrNeighbours ++;
      	neighbours[2] = '1';
      }
      if (  (int)image( Z2i::Point(i,   j+1)) > 0) {
      	nbrNeighbours ++;
      	neighbours[3] = '1';
      }
      if (  (int)image( Z2i::Point(i+1, j+1)) > 0) {
      	nbrNeighbours ++;
      	neighbours[4] = '1';
      }
      if (  (int)image( Z2i::Point(i+1, j  )) > 0) {
      	nbrNeighbours ++;
      	neighbours[5] = '1';
      }
      if (  (int)image( Z2i::Point(i+1, j-1)) > 0) {
      	nbrNeighbours ++;
      	neighbours[6] = '1';
      }
      if (  (int)image( Z2i::Point(i,   j-1)) > 0) {
      	nbrNeighbours ++;
      	neighbours[7] = '1';
      }

      switch (nbrNeighbours) {
      	case 0: // it's an issue: our skeleton is not connex, due to noise, we keep the point
      		break;
      	case 1: // extremity point: we keep it
      		break;
      	case 2: // if accepted we erase it, else we keep it
      		if (isInside(accepted2, 8, neighbours))
      			image.setValue(skeleton[ind], 2);
      		break;
      	case 3:
      		if (isInside(accepted3, 8, neighbours))
      			image.setValue(skeleton[ind], 2);
      		break;
      	case 4:
      		if (isInside(accepted4, 8, neighbours))
      			image.setValue(skeleton[ind], 2);
      		break;
      	case 5: 
      		if (isInside(accepted5, 8, neighbours))
      			image.setValue(skeleton[ind], 2);
      		break;
      	case 6: 
      		if (isInside(accepted6, 8, neighbours))
      			image.setValue(skeleton[ind], 2);
      		break;
      	case 7: // we can erase it
    			image.setValue(skeleton[ind], 2);
      		break;
      	case 8: // inside the shape: we don't do anything 
      		break;
      	default: // can not happen
      		break;
      }
		}

		for (int ind = 0; ind < skeleton.size(); ind++) {
			if ((int)(image(skeleton[ind])) == 2) {
				image.setValue(skeleton[ind], 0);
				skeleton.erase(skeleton.begin()+ind);
				ind--;
				changeMade = true;
			}
		}
	}
}

/*
void ThinSubiteration1(Image2D & pSrc, Image2D & pDst) {
	int rows = pSrc.domain().upperBound()[0];
	int cols = pSrc.domain().upperBound()[1];
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			pDst.setValue(Point(i, j), pSrc(Point(i, j)));

	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			if((int)pSrc(Point(i, j)) == 255) {
				/// get 8 neighbors
				/// calculate C(p)
				int neighbor0 = (int) (pSrc(Point( i-1, j-1)) > 0);
				int neighbor1 = (int) (pSrc(Point( i-1, j)) > 0);
				int neighbor2 = (int) (pSrc(Point( i-1, j+1)) > 0);
				int neighbor3 = (int) (pSrc(Point( i, j+1)) > 0);
				int neighbor4 = (int) (pSrc(Point( i+1, j+1)) > 0);
				int neighbor5 = (int) (pSrc(Point( i+1, j)) > 0);
				int neighbor6 = (int) (pSrc(Point( i+1, j-1)) > 0);
				int neighbor7 = (int) (pSrc(Point( i, j-1)) > 0);
				int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
					 int(~neighbor3 & ( neighbor4 | neighbor5)) +
					 int(~neighbor5 & ( neighbor6 | neighbor7)) +
					 int(~neighbor7 & ( neighbor0 | neighbor1));
				if(C == 1) {
					/// calculate N
					int N1 = int(neighbor0 | neighbor1) +
						int(neighbor2 | neighbor3) +
						int(neighbor4 | neighbor5) +
						int(neighbor6 | neighbor7);
					int N2 = int(neighbor1 | neighbor2) +
						int(neighbor3 | neighbor4) +
						int(neighbor5 | neighbor6) +
						int(neighbor7 | neighbor0);
					int N = min(N1,N2);
					if ((N == 2) || (N == 3)) {
						/// calculate criteria 3
						int c3 = ( neighbor1 | neighbor2 | ~neighbor4) & neighbor3;
						if(c3 == 0) {
							pDst.setValue(Point( i, j), 0);
						}
					}
				}
			}
		}
	}
}


void ThinSubiteration2(Image2D & pSrc, Image2D & pDst) {
	int rows = pSrc.domain().upperBound()[0];
	int cols = pSrc.domain().upperBound()[1];
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			pDst.setValue(Point(i, j), pSrc(Point(i, j)));

	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			if ((int)pSrc(Point( i, j)) > 0) {
				/// get 8 neighbors
				/// calculate C(p)
				int neighbor0 = (int) (pSrc(Point( i-1, j-1)) > 0);
				int neighbor1 = (int) (pSrc(Point( i-1, j)) > 0);
				int neighbor2 = (int) (pSrc(Point( i-1, j+1)) > 0);
				int neighbor3 = (int) (pSrc(Point( i, j+1)) > 0);
				int neighbor4 = (int) (pSrc(Point( i+1, j+1)) > 0);
				int neighbor5 = (int) (pSrc(Point( i+1, j)) > 0);
				int neighbor6 = (int) (pSrc(Point( i+1, j-1)) > 0);
				int neighbor7 = (int) (pSrc(Point( i, j-1)) > 0);
				int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
					int(~neighbor3 & ( neighbor4 | neighbor5)) +
					int(~neighbor5 & ( neighbor6 | neighbor7)) +
					int(~neighbor7 & ( neighbor0 | neighbor1));
				if(C == 1) {
					/// calculate N
					int N1 = int(neighbor0 | neighbor1) +
						int(neighbor2 | neighbor3) +
						int(neighbor4 | neighbor5) +
						int(neighbor6 | neighbor7);
					int N2 = int(neighbor1 | neighbor2) +
						int(neighbor3 | neighbor4) +
						int(neighbor5 | neighbor6) +
						int(neighbor7 | neighbor0);
					int N = min(N1,N2);
					if((N == 2) || (N == 3)) {
						int E = (neighbor5 | neighbor6 | ~neighbor0) & neighbor7;
						if(E == 0) {
							pDst.setValue(Point(i, j), 0);
						}
					}
				}
			}
		}
	}
}


void getSkeleton(Image2D & inputarray, Image2D & outputarray) {
	bool bDone = false;
	int rows = inputarray.domain().upperBound()[0];
	int cols = inputarray.domain().upperBound()[1];


	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			outputarray.setValue(Point(i, j), inputarray(Point(i, j)));



  /// start to thin
	Image2D p_thinMat1 (inputarray);
	Image2D p_thinMat2 (inputarray);

	int turns(0);
  while (bDone != true) {
  	cout << ++turns << endl;
  	bDone = true;
    /// sub-iteration 1
    ThinSubiteration1(inputarray, p_thinMat1);
    /// sub-iteration 2
    ThinSubiteration2(p_thinMat1, p_thinMat2);
    /// compare
    for (int i = 0; i < rows; i++) {
    	for (int j = 0; j < cols; j++) {
    		if ((int)(inputarray(Point(i,j))) != (int)(p_thinMat2(Point(i,j)))) {
    			bDone = false;
    			inputarray.setValue(Point(i,j), (int)(p_thinMat2(Point(i,j))));
    		}
    	}
    }

  }
  // copy result
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      outputarray.setValue(Point( i, j), (int)inputarray(Point( i, j)));
    }
  }
}
*/


void ThinSubiteration1(Image2D & src, vector<Point> &pSrc) {
	int rows = src.domain().upperBound()[0];
	int cols = src.domain().upperBound()[1];

	for (int ind = 0; ind < pSrc.size(); ind++) {
		int i = pSrc[ind][0];
		int j = pSrc[ind][1];
		if((int)src(Point(i, j)) == 255) {
			/// get 8 neighbors
			/// calculate C(p)
			int neighbor0 = (int) (src(Point( i-1, j-1)) > 0);
			int neighbor1 = (int) (src(Point( i-1, j)) > 0);
			int neighbor2 = (int) (src(Point( i-1, j+1)) > 0);
			int neighbor3 = (int) (src(Point( i, j+1)) > 0);
			int neighbor4 = (int) (src(Point( i+1, j+1)) > 0);
			int neighbor5 = (int) (src(Point( i+1, j)) > 0);
			int neighbor6 = (int) (src(Point( i+1, j-1)) > 0);
			int neighbor7 = (int) (src(Point( i, j-1)) > 0);
			int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
				 int(~neighbor3 & ( neighbor4 | neighbor5)) +
				 int(~neighbor5 & ( neighbor6 | neighbor7)) +
				 int(~neighbor7 & ( neighbor0 | neighbor1));
			if(C == 1) {
				/// calculate N
				int N1 = int(neighbor0 | neighbor1) +
					int(neighbor2 | neighbor3) +
					int(neighbor4 | neighbor5) +
					int(neighbor6 | neighbor7);
				int N2 = int(neighbor1 | neighbor2) +
					int(neighbor3 | neighbor4) +
					int(neighbor5 | neighbor6) +
					int(neighbor7 | neighbor0);
				int N = min(N1,N2);
				if ((N == 2) || (N == 3)) {
					/// calculate criteria 3
					int c3 = ( neighbor1 | neighbor2 | ~neighbor4) & neighbor3;
					if(c3 == 0) {
						src.setValue(Point( i, j), 2);
					}
				}
			}
		}
	}
	for (int ind = 0; ind < pSrc.size(); ind++) {
		int i = pSrc[ind][0];
		int j = pSrc[ind][1];
		if((int)src(Point(i, j)) == 2) {
			src.setValue(Point(i, j), 0);
			pSrc.erase(pSrc.begin()+ind);
			ind--;
		}
	}
}


void ThinSubiteration2(Image2D & src, vector<Point> &pSrc) {
	int rows = src.domain().upperBound()[0];
	int cols = src.domain().upperBound()[1];

	for (int ind = 0; ind < pSrc.size(); ind++) {
		int i = pSrc[ind][0];
		int j = pSrc[ind][1];
		if ((int)src(Point( i, j)) > 0) {
			/// get 8 neighbors
			/// calculate C(p)
			int neighbor0 = (int) (src(Point( i-1, j-1)) > 0);
			int neighbor1 = (int) (src(Point( i-1, j)) > 0);
			int neighbor2 = (int) (src(Point( i-1, j+1)) > 0);
			int neighbor3 = (int) (src(Point( i, j+1)) > 0);
			int neighbor4 = (int) (src(Point( i+1, j+1)) > 0);
			int neighbor5 = (int) (src(Point( i+1, j)) > 0);
			int neighbor6 = (int) (src(Point( i+1, j-1)) > 0);
			int neighbor7 = (int) (src(Point( i, j-1)) > 0);
			int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
				int(~neighbor3 & ( neighbor4 | neighbor5)) +
				int(~neighbor5 & ( neighbor6 | neighbor7)) +
				int(~neighbor7 & ( neighbor0 | neighbor1));
			if(C == 1) {
				/// calculate N
				int N1 = int(neighbor0 | neighbor1) +
					int(neighbor2 | neighbor3) +
					int(neighbor4 | neighbor5) +
					int(neighbor6 | neighbor7);
				int N2 = int(neighbor1 | neighbor2) +
					int(neighbor3 | neighbor4) +
					int(neighbor5 | neighbor6) +
					int(neighbor7 | neighbor0);
				int N = min(N1,N2);
				if((N == 2) || (N == 3)) {
					int E = (neighbor5 | neighbor6 | ~neighbor0) & neighbor7;
					if(E == 0) {
						src.setValue(Point( i, j), 2);
					}
				}
			}
		}
	}
	for (int ind = 0; ind < pSrc.size(); ind++) {
		int i = pSrc[ind][0];
		int j = pSrc[ind][1];
		if((int)src(Point(i, j)) == 2) {
			src.setValue(Point(i, j), 0);
			pSrc.erase(pSrc.begin()+ind);
			ind--;
		}
	}
}


void getSkeleton(Image2D & inputarray) {
	bool bDone = false;
	int rows = inputarray.domain().upperBound()[0];
	int cols = inputarray.domain().upperBound()[1];



	vector<Point> skeleton;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			if ((int)inputarray(Point(i,j)) > 0)
				skeleton.push_back(Point(i,j));


  /// start to thin

	int turns(0);
  while (bDone != true) {
  	cout << ++turns << endl;
  	bDone = true;
  	int begin_size (skeleton.size());
    /// sub-iteration 1
    ThinSubiteration1(inputarray, skeleton);
    /// sub-iteration 2
    ThinSubiteration2(inputarray, skeleton);
    /// compare

    if (skeleton.size() != begin_size)
    	bDone = false;
  }
}

void invariants_image(Image2D &image, vector<double> &inv, vector<Point> skeleton) {
	vector<double> temp (rappPerimAire(image));
	inv.push_back(temp[0]);
	inv.push_back(temp[1]);
	getSkeleton(image, skeleton);


	double nbrAnkle (0), nbrExtrems (0), nbrPoints (0);
	for (int ind = 0; ind < skeleton.size(); ind++) {
		int i (skeleton[ind][0]), j (skeleton[ind][1]);
		int nbrNeighbours (0);
    if (  (int)image( Z2i::Point(i-1, j-1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i-1, j  )) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i-1, j+1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i,   j+1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i+1, j+1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i+1, j  )) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i+1, j-1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i,   j-1)) > 0) {
    	nbrNeighbours ++;
    }
    if (nbrNeighbours == 1)
    	nbrExtrems++;
    else if (nbrNeighbours > 2)
    	nbrAnkle++;
	}
	inv.push_back((double)skeleton.size() / (double)image.size());
	inv.push_back(nbrExtrems);
	inv.push_back(nbrAnkle);
}

/*
void invariants_image(Image2D &image, vector<double> &inv, vector<Z2i::Point> skeleton) {
	vector<double> temp (rappPerimAire(image));
	inv.push_back(temp[0]);
	inv.push_back(temp[1]);
	getSkeleton(image, skeleton);
	inv.push_back((double)skeleton.size() / (( image.domain().upperBound()[0] * image.domain().upperBound()[1]) ) );

	double nbrAnkle (0), nbrExtrems (0);
	for (int ind = 0; ind < skeleton.size(); ind++) {
		int i (skeleton[ind][0]), j (skeleton[ind][1]);
		int nbrNeighbours (0);
    if (  (int)image( Z2i::Point(i-1, j-1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i-1, j  )) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i-1, j+1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i,   j+1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i+1, j+1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i+1, j  )) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i+1, j-1)) > 0) {
    	nbrNeighbours ++;
    }
    if (  (int)image( Z2i::Point(i,   j-1)) > 0) {
    	nbrNeighbours ++;
    }
    if (nbrNeighbours == 1)
    	nbrExtrems++;
    else if (nbrNeighbours > 2)
    	nbrAnkle++;
	}
	inv.push_back(nbrExtrems);
	inv.push_back(nbrAnkle);
}*/





//                                                     ____________________________________
//                                                    |                                    |
//                                                    |                 main               |
//                                                    |____________________________________|

int main(int argc, char *argv[])
{

	string filename = (argc < 2) ? "../database/bat-5.pgm" : argv[1];


	//We read the PGM file
	Image2D image = PGMReader<Image2D>::importPGM(filename);
	//trace.info() << "Image2D read :"<< image << endl;
	opening(image);
	closing(image);

	//We convert pixels in ]0,255] into a digital set
	DigitalSet set2d( image.domain() );
	for (int i = 0; i < image.domain().upperBound()[0]; i++)
		for (int j = 0; j < image.domain().upperBound()[1]; j++)
			if ((int)(image(Z2i::Point(i,j))) > 0)
				set2d.insert(Z2i::Point(i,j));

	vector<double> inv = invariants(set2d);
	vector<Point> skeleton;


  invariants_image(image, inv, skeleton);
	cout << inv[0];
	for(int i=1; i<inv.size(); i++)
	{
		cout << "," << inv[i];
	}

	return 0;
}
