#ifndef SHORT_CALCULATIONS_STRUCT_H
#define SHORT_CALCULATIONS_STRUCT_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/determinant.h>

/*
	The structures which calculate the polynomials needed for failure times.
	They immitate the cartesian_predicates_2/3.h (in CGAL/Kinetic/internal/Kernel)
	files, which define them in the same manner - a dummy default constructor, 
	and a single function - operator(). These structures are used in the manner 
	of functions to calculate the function (in most cases it's a polynomial) that
	defines the sign of the short certificates, as time changes.
*/

template <class KK>
struct ShortEdgeCheck2D
{
	ShortEdgeCheck2D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_2 Argument;

	result_type operator()(Argument a, Argument b, double squaredAlpha) const
	{
		typedef typename KK::Certificate_function FT;
		
		FT ax = a.x();
		FT ay = a.y();
		FT bx = b.x();
		FT by = b.y();
        
		//The radius of the smallest circumcircle around the edge is half the length of the edge
		//(a fourth of the squared length, also need to compare it against the squared alpha)
		return squaredAlpha - (CGAL::square(ax - bx) + CGAL::square(ay - by))/4;
	}
};

template<class KK>
struct ShortTriangleCheck2D
{
	ShortTriangleCheck2D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_2 Argument;

	result_type operator()(Argument a, Argument b, Argument c, double squaredAlpha) const
	{

		//The radius of the circumcircle of ABC is |AB||BC||AC| / 2|ABC|
		//where |ABC| is the face of the triangle
		//To avoid the square root, we consider the squared lengths of the sides.

		//Length of the sides (squared)
		Argument distanceAB = CGAL::square(a.x() - b.x()) + CGAL::square(a.y() - b.y());
		Argument distanceBC = CGAL::square(c.x() - b.x()) + CGAL::square(c.y() - b.y());
		Argument distanceAC = CGAL::square(a.x() - c.x()) + CGAL::square(a.y() - c.y());

		/* The face of the triangle in 2D is M^{a,b,c}_[0,1,2]

			|1 ax ay|
 			|1 bx by| = bx*cy + cx*ay + ax*by - (bx*ay + cx*by + ax*cy)
			|1 cx cy|
		*/

		Argument faceABC = b.x() * c.y() + c.x() * a.y() + a.x() * b.y() -
						  (b.x() * a.y() + c.x() * b.y() + a.x() * c.y());

		return squaredAlpha - distanceAB * distanceBC * distanceBC / (4 * CGAL::square(faceABC));
	}
};

///A structure for calculating the certificate functions for short edges in 3D.
///
template <class KK>
struct ShortEdgeCheck3D
{
	ShortEdgeCheck3D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_3 Argument;
    typedef typename KK::Motion_function Motion_function;

	result_type operator()(Argument a, Argument b, Motion_function squaredAlpha) const
	{
		result_type dx = a.x() - b.x();
		result_type dy = a.y() - b.y();
		result_type dz = a.z() - b.z();

        result_type ret = squaredAlpha * 4 - (dx*dx + dy*dy + dz*dz);
		return ret;
	}
};

///A structure for calculating the certificate functions for short faces in 3D.
///
template <class KK>
struct ShortTriangleCheck3D
{
	ShortTriangleCheck3D(){}

	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_3 Argument;
    typedef typename KK::Motion_function Motion_function;

	/// Calculate the 3D polynomial 
	result_type operator()(Argument a, Argument b, Argument c, Motion_function squaredAlpha) const
	{
		
		//The radius of the circumcircle of ABC is |AB||BC||AC| / 2|ABC|
		//where |ABC| is the face of the triangle
		//To avoid computing square root, we consider the squared expression.

		//Length of the sides (squared)
        
        result_type axbx = a.x() - b.x();
        result_type ayby = a.y() - b.y();
        result_type azbz = a.z() - b.z();
        result_type cxbx = c.x() - b.x();
        result_type cyby = c.y() - b.y();
        result_type czbz = c.z() - b.z();
        result_type axcx = a.x() - c.x();
        result_type aycy = a.y() - c.y();
        result_type azcz = a.z() - c.z();

		result_type distanceAB = axbx*axbx + ayby*ayby + azbz*azbz;
		result_type distanceBC = cxbx*cxbx + cyby*cyby + czbz*czbz;
		result_type distanceAC = axcx*axcx + aycy*aycy + azcz*azcz;

		/* 
           The face of the triangle in 3D is 
		   sqrt ((M^{a,b,c}_[0,1,2])^2 + (M^{a,b,c}_[0,2,3])^2 + (M^{a,b,c}_[0,1,3])^2) / 2
		
								|1 ax ay|
 			M^{a,b,c}_[0,1,2] = |1 bx by| = bx*cy + cx*ay + ax*by - (bx*ay + cx*by + ax*cy)
								|1 cx cy|

								|1 ay az|
 			M^{a,b,c}_[0,2,3] = |1 by bz| = by*cz + cy*az + ay*bz - (by*az + cy*bz + ay*cz)
								|1 cy cz|

								|1 ax az|
 			M^{a,b,c}_[0,1,3] = |1 bx bz| = bx*cz + cx*az + ax*bz - (bx*az + cx*bz + ax*cz)
								|1 cx cz|
		*/
		
		result_type M012 = b.x() * c.y() + c.x() * a.y() + a.x() * b.y() -
					      (b.x() * a.y() + c.x() * b.y() + a.x() * c.y());

		result_type M023 = b.y() * c.z() + c.y() * a.z() + a.y() * b.z() -
					      (b.y() * a.z() + c.y() * b.z() + a.y() * c.z());

		result_type M013 = b.x() * c.z() + c.x() * a.z() + a.x() * b.z() -
					      (b.x() * a.z() + c.x() * b.z() + a.x() * c.z());


		return squaredAlpha * (M012 * M012 + M023*M023 + M013*M013) * 4 - 
                    distanceAB * distanceBC * distanceAC;
	}
};

///A structure for calculating the certificate functions for short cells in 3D.
///
template <class KK>
struct ShortTetrahedronCheck3D
{
	ShortTetrahedronCheck3D(){
        std::vector<double> oneV(1, 1);
        one = result_type(oneV.begin(), oneV.end());
    }

	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_3 Argument;
    typedef typename KK::Motion_function Motion_function;

    result_type one;

	/// Calculate the 3D polynomial 
	result_type operator()(Argument a, Argument b, Argument c, Argument d, Motion_function squaredAlpha) const
	{
		result_type SizeA = a.x() * a.x() + a.y() * a.y() + a.z() * a.z();
		result_type SizeB = b.x() * b.x() + b.y() * b.y() + b.z() * b.z();
		result_type SizeC = c.x() * c.x() + c.y() * c.y() + c.z() * c.z();
		result_type SizeD = d.x() * d.x() + d.y() * d.y() + d.z() * d.z();
		
		result_type M1240 =
			CGAL::determinant(a.x(), a.y(), SizeA, one,
							  b.x(), b.y(), SizeB, one,
							  c.x(), c.y(), SizeC, one,
							  d.x(), d.y(), SizeD, one);
		
		result_type M1340 =
			CGAL::determinant(a.x(), a.z(), SizeA, one,
							  b.x(), b.z(), SizeB, one,
							  c.x(), c.z(), SizeC, one,
							  d.x(), d.z(), SizeD, one);
		
		result_type M2340 =
			CGAL::determinant(a.y(), a.z(), SizeA, one,
							  b.y(), b.z(), SizeB, one,
							  c.y(), c.z(), SizeC, one,
							  d.y(), d.z(), SizeD, one);

		result_type M1234 =
			CGAL::determinant(a.x(), a.y(), a.z(), SizeA,
							  b.x(), b.y(), b.z(), SizeB,
							  c.x(), c.y(), c.z(), SizeC,
							  d.x(), d.y(), d.z(), SizeD);

		result_type M1230 = 
			CGAL::determinant(a.x(), a.y(), a.z(), one,
							  b.x(), b.y(), b.z(), one,
							  c.x(), c.y(), c.z(), one,
							  d.x(), d.y(), d.z(), one);
		
		return squaredAlpha * M1230 * M1230 * 4 - 
              (M1240 * M1240 + M1340 * M1340 + M2340 * M2340 - M1234*M1230*4);
	}
};


///Class used for calculating 
template <class KK>
struct GabrielEdgeCheck2D
{
	GabrielEdgeCheck2D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_2 Argument;
	

	///Get the polynomial responsible for the sign change of the Gabriel certificate for edge (a, b)
	result_type operator()(Argument a, Argument b, Argument c) 
	{
		//The center of the smallest circumcircle(the middle of the edge)
		Argument center((a.x() + b.x())/2,(a.y() + b.y())/2);
		
		//The distances between the center and the point that is tested
		Argument dx = c.x() - center.x();
		Argument dy = c.y() - center.y();
		
		//Used for calculating the radius of the circumcircle
		Argument dx1 = a.x() - b.x();
		Argument dy1 = a.y() - b.y();

		//The distance to the point - the radius of the smallest circumcircle around the edge

		return (CGAL::square(dx) + CGAL::square(dy)) - (CGAL::square(dx1) + CGAL::square(dy1))/4;
	}
};

template <class KK>
struct GabrielEdgeCheck3D
{
	GabrielEdgeCheck3D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_3 Argument;
	

	///Get the polynomial responsible for the sign change of the Gabriel certificate for edge (a, b)
	result_type operator()(Argument a, Argument b, Argument c)
	{
		//The center of the smallest circumcircle(the middle of the edge)
		Argument center((a.x() + b.x())/2,(a.y() + b.y())/2,(a.z() + b.z())/2);
		
		//The distances between the center and the point that is tested
		Argument dx = c.x() - center.x();
		Argument dy = c.y() - center.y();
		Argument dz = c.y() - center.z();
		
		//Used for calculating the radius of the circumcircle
		Argument dx1 = a.x() - b.x();
		Argument dy1 = a.y() - b.y();
		Argument dz1 = a.y() - b.z();

		//The distance to the point - the radius of the smallest circumcircle around the edge

		return (CGAL::square(dx) + CGAL::square(dy)+ CGAL::square(dz)) - (CGAL::square(dx1) + CGAL::square(dy1)+ CGAL::square(dz1))/4;
	}
};

template <class KK>
struct GabrielTriangleCheck3D
{
	
	GabrielTriangleCheck3D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_3 Argument;
	
	result_type operator()(Argument a, Argument b, Argument c, Argument d)
	{
		Argument SizeA = CGAL::square(a.x()) + CGAL::square(a.y()) + CGAL::square(a.z());
		Argument SizeB = CGAL::square(b.x()) + CGAL::square(b.y()) + CGAL::square(b.z());
		Argument SizeC = CGAL::square(c.x()) + CGAL::square(c.y()) + CGAL::square(c.z());
		Argument SizeD = CGAL::square(d.x()) + CGAL::square(d.y()) + CGAL::square(d.z());

		result_type M2340 =
			CGAL::determinant(a.y(), a.z(), SizeA, 1,
							  b.y(), b.z(), SizeB, 1,
							  c.y(), c.z(), SizeC, 1,
							  d.y(), d.z(), SizeD, 1);

		result_type M230 =
			CGAL::determinant(a.y(), a.z(), 1,
							  b.y(), b.z(), 1,
							  c.y(), c.z(), 1);

		result_type M1340 =
			CGAL::determinant(a.x(), a.z(), SizeA, 1,
							  b.x(), b.z(), SizeB, 1,
							  c.x(), c.z(), SizeC, 1,
							  d.x(), d.z(), SizeD, 1);

		result_type M130 =
			CGAL::determinant(a.x(), a.z(), 1,
							  b.x(), b.z(), 1,
							  c.x(), c.z(), 1);

		result_type M1240 =
			CGAL::determinant(a.x(), a.y(), SizeA, 1,
							  b.x(), b.y(), SizeB, 1,
							  c.x(), c.y(), SizeC, 1,
							  d.x(), d.y(), SizeD, 1);

		result_type M120 =
			CGAL::determinant(a.x(), a.y(), 1,
							  b.x(), b.y(), 1,
							  c.x(), c.y(), 1);

		result_type M1230 =
			CGAL::determinant(a.x(), a.y(), a.z(), 1,
							  b.x(), b.y(), b.z(), 1,
							  c.x(), c.y(), c.z(), 1,
							  d.x(), d.y(), d.z(), 1);

		result_type M123 =
			CGAL::determinant(a.x(), a.y(), a.z(),
							  b.x(), b.y(), b.z(),
							  c.x(), c.y(), c.z());

		return (M2340 * M230) + (M1340 * M130) + (M1240 * M120) - 2*(M1230 * M123)
	}

};

#endif