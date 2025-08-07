/*! \file lineseg.h
 
\brief Classes for representing a line segment
\author Jyh-Ming Lien
\author <a href="http://masc.cs.gmu.edu">MASC group</a>, George Mason University
*/

#pragma once
#ifndef _CUSG_LIENSEG_H_
#define _CUSG_LIENSEG_H_

#include "Point.h"

namespace cusg
{
    /// a line segment
    //template<class T>
    class LineSeg2d
    {
    public:

        LineSeg2d(){}
        LineSeg2d(const LineSeg2d& other){ pt[0]=other.pt[0]; pt[1]=other.pt[1]; }
        LineSeg2d(double p1x, double p1y, double p2x, double p2y){ pt[0].set(p1x,p1y); pt[1].set(p2x,p2y); }
        LineSeg2d(const Point2d& t1, const Point2d& t2){ pt[0]=t1; pt[1]=t2; }

        Point2d pt[2]; //!< end points of the line segment

    };

    class GIS_LineSeg2d : public LineSeg2d
    {
    public:

        GIS_LineSeg2d() : LineSeg2d()
        { m_seg_id=m_curve_id=0; m_start_geo_dist=m_end_geo_dist=0; }

        GIS_LineSeg2d(const GIS_LineSeg2d& other): LineSeg2d(other)
        {
            m_seg_id=other.m_seg_id; m_curve_id=other.m_curve_id;
            m_end_geo_dist=other.m_end_geo_dist;
            m_start_geo_dist=other.m_start_geo_dist;
        }

        GIS_LineSeg2d(double p1x, double p1y, double p2x, double p2y) : LineSeg2d(p1x,p1y,p2x,p2y)
        { m_seg_id=m_curve_id=0; m_start_geo_dist=m_end_geo_dist=0;}

        GIS_LineSeg2d(const Point2d& t1, const Point2d& t2):LineSeg2d(t1,t2)
        { m_seg_id=m_curve_id=0; m_start_geo_dist=m_end_geo_dist=0; }

        double m_start_geo_dist;
        double m_end_geo_dist;
        int m_seg_id;   //id of this segment
        int m_curve_id; //id of the curve (polyline) that this segment belongs

		void rotate(double cx, double cy, double angle)
		{
			pt[0].rotate(cx, cy, angle);
			pt[1].rotate(cx, cy, angle);
		}

    };

}//namespace cusg

#endif //_CUSG_LIENSEG_H_


