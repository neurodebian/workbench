/*LICENSE_START*/
/*
 *  Copyright (C) 2014  Washington University School of Medicine
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*LICENSE_END*/

#include <cmath>

#include "CaretLogger.h"
#include "MathFunctions.h"
#include "Plane.h"
#include "Vector3D.h"

using namespace caret;

/**
 * Static method
 * @return a plane that is in a string created by Plane::toFormattedString()
 * @ s
 *    The string
 * @sealso toFormattedString
 */
Plane
Plane::fromFormattedString(const AString& s)
{
    std::vector<float> v;
    s.toNumbers(s, v);
    
    if (v.size() == 8) {
        if (v[7] != 0.0) {
            const Vector3D pointOnPlane(v[4], v[5], v[6]);
            return Plane(v[0], v[1], v[2], v[3],
                         pointOnPlane);
        }
    }
    
    CaretLogSevere("ERROR: Unable to decode Plane from string \""
                   + s
                   + "\"");
    return Plane();
}

/**
 * Construct an invalid plane.
 * Intended for use by the assignement operator.
 */
Plane::Plane()
: CaretObject()
{
    m_pointOnPlane[0] = 0.0;
    m_pointOnPlane[1] = 0.0;
    m_pointOnPlane[2] = 0.0;
    
    m_normalVector[0] = 0.0;
    m_normalVector[1] = 0.0;
    m_normalVector[2] = 0.0;

    
    m_validPlaneFlag = false;
    
    m_A = m_normalVector[0];
    m_B = m_normalVector[1];
    m_C = m_normalVector[2];
    m_D = 0.0;
}

/**
 * Construct a plane from three points that are on the plane.
 * These points should be in counter-clockwise order.
 *
 * @param p1 Point on plane.
 * @param p2 Point on plane.
 * @param p3 Point on plane.
 *
 */
Plane::Plane(const float p1[3],
             const float p2[3],
             const float p3[3])
    : CaretObject()
{
    m_pointOnPlane[0] = p1[0];
    m_pointOnPlane[1] = p1[1];
    m_pointOnPlane[2] = p1[2];
    
    double p2d[3] = {
        p2[0],
        p2[1],
        p2[2]
    };
    
    double p3d[3] = {
        p3[0],
        p3[1],
        p3[2]
    };
    
    m_validPlaneFlag =  MathFunctions::normalVector(m_pointOnPlane,
                                                        p2d,
                                                        p3d,
                                                        m_normalVector);
    
    //
    // Compute the plane equation
    //
    m_A = m_normalVector[0];
    m_B = m_normalVector[1];
    m_C = m_normalVector[2];
    m_D = -(m_A*p1[0] + m_B*p1[1] + m_C*p1[2]);
}

/**
 * Construct a plane from a unit normal vector (length = 1) and a point on the plane.
 *
 * @param normalVector
 *     The normal vector of the plane.
 * @param pointOnPlane
 *     A point on the plane.
 */
Plane::Plane(const float unitNormalVector[3],
             const float pointOnPlane[3])
{
    m_pointOnPlane[0] = pointOnPlane[0];
    m_pointOnPlane[1] = pointOnPlane[1];
    m_pointOnPlane[2] = pointOnPlane[2];
    
    m_normalVector[0] = unitNormalVector[0];
    m_normalVector[1] = unitNormalVector[1];
    m_normalVector[2] = unitNormalVector[2];

    m_A = m_normalVector[0];
    m_B = m_normalVector[1];
    m_C = m_normalVector[2];    
    m_D = (-m_A * m_pointOnPlane[0]
               -m_B * m_pointOnPlane[1]
               -m_C * m_pointOnPlane[2]);

    m_validPlaneFlag = (MathFunctions::vectorLength(m_normalVector) > 0.0);
}

/**
 * Contructor from plane equation and a point on the plane
 * @param A
 *   A component of plane equation
 * @param B
 *   B component of plane equation
 * @param C
 *   C component of plane equation
 * @param D
 *   D component of plane equation
 * @param pointOnPlane
 *   Point on the plane
 */
Plane::Plane(const float A,
             const float B,
             const float C,
             const float D,
             const float pointOnPlane[3])
{
    m_pointOnPlane[0] = pointOnPlane[0];
    m_pointOnPlane[1] = pointOnPlane[1];
    m_pointOnPlane[2] = pointOnPlane[2];
    
    m_normalVector[0] = A;
    m_normalVector[1] = B;
    m_normalVector[2] = C;
    
    m_A = A;
    m_B = B;
    m_C = C;
    m_D = D;

    m_validPlaneFlag = (MathFunctions::vectorLength(m_normalVector) > 0.0);
}

/**
 * Copy constructor.
 * @param p
 *    Object that is copied.
 */
Plane::Plane(const Plane& p)
: CaretObject(p)
{
    this->copyHelperPlane(p);
}

/**
 * Assignment operator.
 * @param p
 *    Data copied from obj to this.
 * @return
 *    Reference to this object.
 */
Plane&
Plane::operator=(const Plane& p)
{
    if (this != &p) {
        CaretObject::operator=(p);
        this->copyHelperPlane(p);
    }
    return *this;
}

/**
 * Helps with copying an object of this type.
 * @param p
 *    Object that is copied.
 */
void
Plane::copyHelperPlane(const Plane& p)
{
    m_A = p.m_A;
    m_B = p.m_B;
    m_C = p.m_C;
    m_D = p.m_D;
    
    m_normalVector[0] = p.m_normalVector[0];
    m_normalVector[1] = p.m_normalVector[1];
    m_normalVector[2] = p.m_normalVector[2];
    
    m_pointOnPlane[0] = p.m_pointOnPlane[0];
    m_pointOnPlane[1] = p.m_pointOnPlane[1];
    m_pointOnPlane[2] = p.m_pointOnPlane[2];
    
    m_validPlaneFlag = p.m_validPlaneFlag;
}

/**
 * Destructor
 */
Plane::~Plane()
{
}

/**
 * Is the plane valid?
 * @return  true if plane is valid, else false.
 *
 */
bool
Plane::isValidPlane() const
{
    return m_validPlaneFlag;
}

/**
 * Find the points where a triangle intersects the plane.
 * @param t1  First point in triangle.
 * @param t2  Second point in triangle.
 * @param t3  Third point in triangle.
 * @param intersectionPointOut1 If there is intersection,
 * this will be loaded with one intersection point.
 * @param intersectionPointOut2 If there is intersection,
 * this will be loaded with other intersection point.
 * @return  true if any two edges of the triangle intersect
 * the plane, else false.
 *
 */
bool
Plane::triangleIntersectPlane(
                   const float t1[3],
                   const float t2[3],
                   const float t3[3],
                   float intersectionPointOut1[3],
                   float intersectionPointOut2[3]) const
{
    float* intersection = intersectionPointOut1;
    int count = 0;
    if (lineSegmentIntersectPlane(t1, t2, intersection)) {
        count++;
        intersection = intersectionPointOut2;
    }
    if (lineSegmentIntersectPlane(t2, t3, intersection)) {
        count++;
        if (count == 2) {
            return true;
        }
        intersection = intersectionPointOut2;
    }
    if (lineSegmentIntersectPlane(t3, t1, intersection)) {
        count++;
        if (count == 2) {
            return true;
        }
    }
    
    return false;
}

/**
 * Get the plane.
 *
 * @param aOut
 *   The value of 'A'
 * @param bOut
 *   The value of 'B'
 * @param cOut
 *   The value of 'C'
 * @param dOut
 *   The value of 'D'
 */
void
Plane::getPlane(double& aOut,
                double& bOut,
                double& cOut,
                double& dOut) const
{
    aOut = m_A;
    bOut = m_B;
    cOut = m_C;
    dOut = m_D;
}

/**
 * @return The plane equation (A, B, C, D) in a 4 element std::array
 */
std::array<double, 4>
Plane::getPlaneEquation() const
{
    std::array<double, 4> eq {
        static_cast<float>(m_A),
        static_cast<float>(m_B),
        static_cast<float>(m_C),
        static_cast<float>(m_D)
    };
    return eq;
}

/**
 * Get the plane's normal vector.
 *
 * @param normalVectorOut
 *   On exit, contains the plane's normal vector.
 */
void
Plane::getNormalVector(double normalVectorOut[3]) const
{
    normalVectorOut[0] = m_normalVector[0];
    normalVectorOut[1] = m_normalVector[1];
    normalVectorOut[2] = m_normalVector[2];
}

/**
 * Get the plane's normal vector.
 *
 * @param normalVectorOut
 *   On exit, contains the plane's normal vector.
 */
void
Plane::getNormalVector(float normalVectorOut[3]) const
{
    normalVectorOut[0] = m_normalVector[0];
    normalVectorOut[1] = m_normalVector[1];
    normalVectorOut[2] = m_normalVector[2];
}

/**
 * @return The plane's normal vector.
 */
Vector3D
Plane::getNormalVector() const
{
    float n[3];
    getNormalVector(n);
    return Vector3D(n);
}

/**
 * Invert the normal vector (multiple by -1)
 */
void
Plane::invertNormalVector()
{
    m_A *= -1.0;
    m_B *= -1.0;
    m_C *= -1.0;
    m_normalVector[0] = m_A;
    m_normalVector[1] = m_B;
    m_normalVector[2] = m_C;
}


/**
 * Get absolute distance of point from the plane.
 * @param p Point.
 * @return  Absolute distance of "p" from plane.
 *
 */
double
Plane::absoluteDistanceToPlane(const float p[3]) const
{
    double dist = (m_A * p[0] + m_B * p[1] + m_C * p[2] + m_D);
    if (dist < 0.0f) dist = -dist;
    return dist;
}

/**
 * Get signed distance of point from the plane.
 * @param p Point.
 * @return  Signed distance of "p" from plane.
 *
 */
double
Plane::signedDistanceToPlane(const float p[3]) const
{
    double dist = (m_A * p[0] + m_B * p[1] + m_C * p[2] + m_D);
    return dist;
}

/**
 * Find the intersection of a line segment with the plane.
 * @param lp1 start of line.
 * @param lp2 end of line.
 * @param intersectionOut If the line intersects the plane,
 * the point of intersection will be loaded into this.  If 
 * there is no intersection, the values of this array will 
 * not be altered.
 * @return  true if the line segment intersects the plane,
 * else false.
 *
 */
bool
Plane::lineSegmentIntersectPlane(
                   const float lp1[3],
                   const float lp2[3],
                   float intersectionOut[3]) const
{
    /*
     * Ray formed by lp1 ==> lp2 (NOT a unit vector; do not normalize)
     */
    double r0 = lp2[0] - lp1[0];
    double r1 = lp2[1] - lp1[1];
    double r2 = lp2[2] - lp1[2];
    
    /*
     * Denominator is dot product of the plane's normal
     * and the ray.  If it is zero, that means the ray
     * is parallel to the plane.
     * 
     * Perhaps should test for value near zero, and,
     * if it is, see if line end points are very 
     * close to plane, and, if so, use the line
     * endpoints, projected to the plane, as the
     * intersection.
     */
    double denominator = 
    (m_A * r0)
    + (m_B * r1) 
    + (m_C * r2);
    if (denominator == 0.0) {
        return false;
    }
    
    /*
     * "t" is the normalized distance of plane from the line origin to
     * the endpoint of the line.  
     * 
     * (0 <= t <= 1): The line intersects the plane with one point
     * above the plane and one point below the plane.
     * 
     * (t > 1): Vector (P1, P2) point towards plane but both are
     * on same side of the plane.
     * 
     * (t < 0): Vector (P1, P2) point away from plane but both are
     * on same side of the plane. 
     */
    double numerator = -((m_A * lp1[0]) + (m_B * lp1[1]) + (m_C * lp1[2]) + m_D);
    double t = numerator / denominator;
    
    if ((t >= 0.0f) && (t <= 1.0f)) {
        intersectionOut[0] = (float)(lp1[0] + r0 * t);
        intersectionOut[1] = (float)(lp1[1] + r1 * t);
        intersectionOut[2] = (float)(lp1[2] + r2 * t);
        return true;
    }
    
    return false;
}

/**
 * Project the given point to the plane.
 * @param pointIn
 *   Point that will be projected.
 * @param pointProjectedOut
 *   Coordinates of point after projection to the plane.
 */
void
Plane::projectPointToPlane(const float pointIn[3],
                           float pointProjectedOut[3]) const
{
    double xo[3] = {
        pointIn[0] - m_pointOnPlane[0],
        pointIn[1] - m_pointOnPlane[1],
        pointIn[2] - m_pointOnPlane[2]
    };
    
    float t = MathFunctions::dotProduct(m_normalVector, xo);
    
    pointProjectedOut[0] = pointIn[0] - (t * m_normalVector[0]);
    pointProjectedOut[1] = pointIn[1] - (t * m_normalVector[1]);
    pointProjectedOut[2] = pointIn[2] - (t * m_normalVector[2]);
}

/**
 * Project the given point to the plane.
 * @param pointIn
 *   Point that will be projected.
 * @return
 *   Coordinates of point after projection to the plane.
 */
Vector3D
Plane::projectPointToPlane(const float pointIn[3]) const
{
    Vector3D pOut;
    projectPointToPlane(pointIn, pOut);
    return pOut;
}

/**
 * Determine if and where a ray intersects the plane.
 *
 * @param rayOrigin
 *     Origin of the ray
 * @param rayVector
 *     Vector defining the ray.
 * @param intersectionOutXYZ
 *     Coordinate of where the ray intersects the plane (XYZ)
 * @param distanceOut
 *     Distance of the ray origin from the plane.
 * @return
 *     True if the ray intersects the plane, else false.
 */
bool
Plane::rayIntersection(const float rayOrigin[3],
                       const float rayVector[3],
                       Vector3D& intersectionOutXYZ,
                       float& distanceOut) const
{
    /* Convert the ray into a unit vector
     *
     */
    double ray[3] = { rayVector[0], rayVector[1], rayVector[2] };
    MathFunctions::normalizeVector(ray);
    
    /*
     * Parametric coordinate of where ray intersects plane
     */
    double denom = m_A * ray[0] + m_B * ray[1] + m_C * ray[2];
    if (denom != 0) {
        const double t = -(m_A * rayOrigin[0] + m_B * rayOrigin[1] + m_C * rayOrigin[2] + m_D) / denom;
        
        intersectionOutXYZ[0] = (float)(rayOrigin[0] + ray[0] * t);
        intersectionOutXYZ[1] = (float)(rayOrigin[1] + ray[1] * t);
        intersectionOutXYZ[2] = (float)(rayOrigin[2] + ray[2] * t);
        
        distanceOut = (float)t;
        
        return true;
    }
    
    return false;
}

/**
 * Shift the plane by the given distance
 * @param distance
 *    Distance to shift the plane
 */
void
Plane::shiftPlane(const float distance)
{
    m_D += distance;
}

/**
 * @return String describing the plane.
 *
 */
AString
Plane::toString() const
{
    AString s;
    if (isValidPlane()) {
        s = (AString::number(m_A)
             + "x + "
             + AString::number(m_B)
             + "y + "
             + AString::number(m_C)
             + "z + "
             + AString::number(m_D));
    }
    else {
        s = "invalid";
    }
    
    return s;
}

/**
 * @return Formatted string for saving
 * @seealso fromFormattedString
 */
AString
Plane::toFormattedString() const
{
    if ( ! isValidPlane()) {
        return "";
    }
    
    /*
     * Default vector to all zeros (invalid plane)
     */
    std::vector<float> v(8, 0.0);
    
    CaretAssertVectorIndex(v, 7);
    v[0] = m_A;
    v[1] = m_B;
    v[2] = m_C;
    v[3] = m_D;
    v[4] = m_pointOnPlane[0];
    v[5] = m_pointOnPlane[1];
    v[6] = m_pointOnPlane[2];
    v[7] = (m_validPlaneFlag ? 1.0 : 0.0);
    
    const AString separator(" ");
    const char    floatingPointFormat('f');
    const int32_t precision(8);
    const AString s(AString::fromNumbers(v,
                                         separator,
                                         floatingPointFormat,
                                         precision));
    
    return s;
}

/**
 * Get text representations of ABCD and point on plane
 * @param abcdOut
 *    Output with ABCD
 * @param pointXyzOut
 *    Output with point on plane
 * @return True if valid, else false.
 */
bool
Plane::toAbcdAndPointXYZ(AString& abcdOut,
                         AString& pointXyzOut) const
{
    if (isValidPlane()) {
        abcdOut = ("A="
                      + AString::number(m_A)
                      + ", B="
                      + AString::number(m_B)
                      +", C="
                      + AString::number(m_C)
                      +", D="
                      + AString::number(m_D));
        pointXyzOut = ("Point on Plane: "
                       + AString::fromNumbers(m_pointOnPlane, 3));
        return true;
    }
    
    abcdOut = "";
    pointXyzOut   = "";
    
    return false;
}


/**
 * Unit test the class.
 *
 */
void
Plane::unitTest(std::ostream& stream,
                const bool /*isVerbose*/)
{
    stream << "Plane::unitTest is starting" << std::endl;
    Plane::unitTest1(stream);
    Plane::unitTest2(stream);
    Plane::unitTest3(stream);
    stream << "Plane::unitTest has ended" << std::endl;
}

/**
 * Unit test a plane and line intersection.
 * @param testName  Name of test.
 * @param p1  Plane Point 1.
 * @param p2  Plane Point 2.
 * @param p3  Plane Point 3.
 * @param l1  Line End Point 1.
 * @param l2  Line End Point 2.
 * @param correctIntersectionPoint
 * @param intersectionValid
 *
 */
void
Plane::unitTestLineIntersectPlane(std::ostream& stream,
                   const AString& testName,
                   const float p1[3],
                   const float p2[3],
                   const float p3[3],
                   const float l1[3],
                   const float l2[3],
                   const float correctIntersectionPoint[3],
                   const bool intersectionValid)
{
    float intersection[3];
    Plane p(p1, p2, p3);
    bool result = p.lineSegmentIntersectPlane(l1, l2, intersection);
    
    AString sb;
    if (intersectionValid && result) {
        if (MathFunctions::compareArrays(correctIntersectionPoint, intersection, 3, 0.001f)) {
            return;
        }
        sb.append("Intersection should be " 
                  + AString::fromNumbers(correctIntersectionPoint, 3, ","));
        sb.append(" but is " + AString::fromNumbers(intersection, 3, ","));
    }
    else if (intersectionValid != result) {
        sb.append("intersection " + AString::fromBool(result) + " but should be " + AString::fromBool(intersectionValid));
    }
    if (sb.length() > 0) {
        sb = ("Line/Plane Intersection Test " + testName + " FAILED: "
              + sb);
        stream << sb.toStdString() << std::endl;
    }
}

/**
 * Are the plane normal vectors orthogonal?
 * @param normalVectorOne
 *    Normal vector from first plane
 * @param normalVectorTwo
 *    Normal vector from second plane
 * @param normalVectorThree
 *    Normal vector from third plane
 * @param optionalMessageOut
 *    Optional (if not NULL) message describing non-orthogonal
 * @return
 *   True if planes appear to be orthogonal, else false
 */
bool
Plane::arePlanesOrthogonal(const Vector3D& normalVectorOne,
                           const Vector3D& normalVectorTwo,
                           const Vector3D& normalVectorThree,
                           AString* optionalMessageOut)
{
    AString errorMessageOut;
    
    /*
     * Note dot product for 90 degrees is -1 or 1
     */
    const float tolerance(0.001);
    const float tolPosZero(tolerance);
    const float tolNegZero(- tolerance);
    const float tolPosOne(1.0 - tolerance);
    const float tolNegOne(-1.0 + tolerance);
    
    /*
     * Dot poduct of zero indicates 90 degree angle between vectors
     */
    const float dotOneTwo(normalVectorOne.dot(normalVectorTwo));
    const float dotOneThree(normalVectorOne.dot(normalVectorThree));
    const float dotTwoThree(normalVectorTwo.dot(normalVectorThree));
    if ((dotOneTwo < tolNegZero)
        || (dotOneTwo > tolPosZero)) {
        const float angleDegrees(MathFunctions::toDegrees(std::acos(MathFunctions::limitRange(dotOneTwo, -1.0f, 1.0f))));
        errorMessageOut.appendWithNewLine("Planes one and two are not 90 degrees; dot product="
                                          + AString::number(dotOneTwo)
                                          + ", angle="
                                          + AString::number(angleDegrees)
                                          + " degrees");
    }

    if ((dotOneThree < tolNegZero)
        || (dotOneThree > tolPosZero)) {
        const float angleDegrees(MathFunctions::toDegrees(std::acos(MathFunctions::limitRange(dotOneThree, -1.0f, 1.0f))));
        errorMessageOut.appendWithNewLine("Planes one and three are not 90 degrees; dot product="
                                          + AString::number(dotOneThree)
                                          + ", angle="
                                          + AString::number(angleDegrees)
                                          + " degrees");
    }

    if ((dotTwoThree < tolNegZero)
        || (dotTwoThree > tolPosZero)) {
        const float angleDegrees(MathFunctions::toDegrees(std::acos(MathFunctions::limitRange(dotTwoThree, -1.0f, 1.0f))));
        errorMessageOut.appendWithNewLine("Planes two and three are not 90 degrees; dot product="
                                          + AString::number(dotTwoThree)
                                          + ", angle="
                                          + AString::number(angleDegrees)
                                          + " degrees");
    }
    
    if (errorMessageOut.isEmpty()) {
        const Vector3D vecOneTwo(normalVectorOne.cross(normalVectorTwo));
        
        const float dot(vecOneTwo.dot(normalVectorThree));
        if ((dot > tolPosOne)
            || (dot < tolNegOne)) {
            return true;
        }
        else {
            const float angleDegrees(MathFunctions::toDegrees(std::acos(MathFunctions::limitRange(dot, -1.0f, 1.0f))));
            errorMessageOut.appendWithNewLine("Vector formed by cross product of first two vectors not aligned with third plane, "
                                              "should be 0 or 180; dot product="
                                              + AString::number(dot)
                                              +", angle="
                                              + AString::number(angleDegrees)
                                              + " degrees.");
        }
    }

    if (optionalMessageOut != NULL) {
        *optionalMessageOut = errorMessageOut;
    }
    return false;
}

/**
 * @return Angle in degrees formed by the normal vectors of the two planes
 * @param p1
 *    First plane
 * @param p2
 *    Second plane
 * @param optionalValidFlagOut
 *    If not NULL, will indicate returned value is valid
 */
float
Plane::angleDegreesOfPlaneNormalVectors(const Plane& p1,
                                        const Plane& p2,
                                        bool* optionalValidFlagOut)
{
    if (optionalValidFlagOut != NULL) {
        *optionalValidFlagOut = false;
    }

    bool validPlanesFlag(true);
    if ( ! p1.isValidPlane()) {
        if (optionalValidFlagOut == NULL) {
            CaretLogSevere("Program error: Plane p1 is invalid.");
        }
        validPlanesFlag = false;
    }
    if ( ! p2.isValidPlane()) {
        if (optionalValidFlagOut == NULL) {
            CaretLogSevere("Program error: Plane p2 is invalid.");
        }
        validPlanesFlag = false;
    }
    
    if ( ! validPlanesFlag) {
        return 0.0;
    }
    
    const Vector3D n1(p1.getNormalVector());
    const Vector3D n2(p2.getNormalVector());
    float dot(n1.dot(n2));
    if (dot <= -1.0) {
        dot = -1.0;
    }
    else if (dot >= 1.0) {
        dot = 1.0;
    }
    const float angle(MathFunctions::toDegrees(std::acos(dot)));
    if (optionalValidFlagOut != NULL) {
        *optionalValidFlagOut = true;
    }
    return angle;
}

void
Plane::unitTest1(std::ostream& stream)
{
    float p1[] = { -50.0f, -60.0f, 2.0f };
    float p2[] = {  50.0f, -60.0f, 2.0f };
    float p3[] = {  50.0f,  60.0f, 2.0f };
    float l1[] = { 25.0f, 10.0f, -3.0f };
    float l2[] = { 25.0f, 10.0f,  7.0f };
    float correctIntersectionPoint[] = {
        25.0f, 10.0f, 2.0f      
    };
    Plane::unitTestLineIntersectPlane(stream, "1", p1, p2, p3, l1, l2, correctIntersectionPoint, true);
}

void
Plane::unitTest2(std::ostream& stream)
{
    float p1[] = { -50.0f, -60.0f, 2.0f };
    float p2[] = {  50.0f, -60.0f, 2.0f };
    float p3[] = {  50.0f,  60.0f, 2.0f };
    float l1[] = { 25.0f, 10.0f,  7.0f };
    float l2[] = { 25.0f, 10.0f,  9.0f };
    float correctIntersectionPoint[3];
    Plane::unitTestLineIntersectPlane(stream, "2", p1, p2, p3, l1, l2, correctIntersectionPoint, false);}

void
Plane::unitTest3(std::ostream& stream)
{
    float p1[] = { -50.0f, -60.0f, 2.0f };
    float p2[] = {  50.0f, -60.0f, 2.0f };
    float p3[] = {  50.0f,  60.0f, 2.0f };
    float l1[] = { 25.0f, 10.0f, -5.0f };
    float l2[] = { 25.0f, 10.0f, -3.0f };
    float correctIntersectionPoint[3];
    Plane::unitTestLineIntersectPlane(stream, "3", p1, p2, p3, l1, l2, correctIntersectionPoint, false);
}

