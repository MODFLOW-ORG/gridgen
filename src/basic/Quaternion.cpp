
#include "Quaternion.h"

namespace cusg{

    //
    //from Jonathan Blow
    //Understanding Slerp, Then Not Using It
    //
    Quaternion slerp(const Quaternion& q0, const Quaternion& q1, REAL t)
    {
        // v0 and v1 should be unit length or else
        // something broken will happen.

        // Compute the cosine of the angle between the two vectors.
        const REAL dot = q0.getReal()*q1.getReal()+q0.getComplex()*q1.getComplex();
        if(dot>=1) return q0;//q0==q1

        const REAL DOT_THRESHOLD = 0.9999995f;
        if (dot > DOT_THRESHOLD) {
            // If the inputs are too close for comfort, linearly interpolate
            // and normalize the result.
            Quaternion result=q0+(q1-q0)*t;
            result.normalize();
            return result;
        }

        //Clamp(dot, -1, 1);                 // Robustness: Stay within domain of acos()
        double theta_0 = acos(dot);          // theta_0 = angle between input vectors
        REAL theta = (REAL)(theta_0*t);    // theta = angle between v0 and result 

        Quaternion q2 = q1-q0*dot;
        q2.normalize();              // { v0, v2 } is now an orthonormal basis

        return q0*cos(theta)+q2*sin(theta);
    }

    istream& operator>>(istream & in, Quaternion & q )
    {
        in>>q.m_s>>q.m_v;
        return in;
    }

    ostream& operator<<(ostream & out, const Quaternion & q )
    {
        out<<q.m_s<<" "<<q.m_v;
        return out;
    }

    Quaternion operator*(const Vector3d & v, const Quaternion & q2)
    {
        Quaternion q1(0,v);
        return q1*q2;
    }
}
