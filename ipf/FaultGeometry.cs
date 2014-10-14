using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using edu.mines.jtk.util;

namespace ipf
{
    public static class FaultGeometry
    {
        /**
        * Returns fault dip vector for specified strike and dip angles.
        * @param phi fault strike angle, in degrees.
        * @param theta fault dip angle, in degrees.
        * @return array {u1,u2,u3} of components for dip vector.
        */
        public static float[] faultDipVectorFromStrikeAndDip(
        double phi, double theta)
        {
            double p = ArrayMath.toRadians(phi);
            double t = ArrayMath.toRadians(theta);
            double cp = Math.Cos(p);
            double sp = Math.Sin(p);
            double ct = Math.Cos(t);
            double st = Math.Sin(t);
            float u1 = (float)(st);
            float u2 = (float)(ct * cp);
            float u3 = (float)(-ct * sp);
            return new float[] { u1, u2, u3 };
        }

        /**
        * Returns fault strike vector for specified strike and dip angles.
        * The dip angle theta is not used, but is provided for consistency.
        * @param phi fault strike angle, in degrees.
        * @param theta fault dip angle, in degrees.
        * @return array {v1,v2,v3} of components for strike vector.
        */
        public static float[] faultStrikeVectorFromStrikeAndDip(
        double phi, double theta)
        {
            double p = ArrayMath.toRadians(phi);
            double cp = Math.Cos(p);
            double sp = Math.Sin(p);
            float v1 = 0.0f;
            float v2 = (float)sp;
            float v3 = (float)cp;
            return new float[] { v1, v2, v3 };
        }

        /**
        * Returns fault normal vector for specified strike and dip angles.
        * @param phi fault strike angle, in degrees.
        * @param theta fault dip angle, in degrees.
        * @return array {w1,w2,w3} of components for normal vector.
        */
        public static float[] faultNormalVectorFromStrikeAndDip(
        double phi, double theta)
        {
            double p = ArrayMath.toRadians(phi);
            double t = ArrayMath.toRadians(theta);
            double cp = Math.Cos(p);
            double sp = Math.Sin(p);
            double ct = Math.Cos(t);
            double st = Math.Sin(t);
            float w1 = (float)(-ct);
            float w2 = (float)(st * cp);
            float w3 = (float)(-st * sp);
            return new float[] { w1, w2, w3 };
        }

        /**
        * Returns fault strike angle for specified fault dip vector.
        * The components u2 and u3 must not both be zero; that is, the
        * fault dip vector cannot be vertical.
        * @param u1 1st component of fault dip vector.
        * @param u2 2nd component of fault dip vector.
        * @param u3 3rd component of fault dip vector.
        * @return fault strike angle, in degrees.
        */
        public static float faultStrikeFromDipVector(float u1, float u2, float u3)
        {
            //Check.argument(u2 != 0.0f || u3 != 0.0f, "dip vector is not vertical");
            return range360(ArrayMath.toDegrees(Math.Atan2(-u3, u2)));
        }

        /**
        * Returns fault strike angle for specified fault dip vector.
        * The components u2 and u3 must not both be zero; that is, the
        * fault dip vector cannot be vertical.
        * @param u array {u1,u2,u3} of components of fault dip vector.
        * @return fault strike angle, in degrees.
        */
        public static float faultStrikeFromDipVector(float[] u)
        {
            return faultStrikeFromDipVector(u[0], u[1], u[2]);
        }

        /**
        * Returns fault dip angle for specified fault dip vector.
        * @param u1 1st component of fault dip vector.
        * @param u2 2nd component of fault dip vector.
        * @param u3 3rd component of fault dip vector.
        * @return fault dip angle, in degrees.
        */
        public static float faultDipFromDipVector(float u1, float u2, float u3)
        {
            return ArrayMath.toDegrees((float)Math.Asin(u1));
        }

        /**
        * Returns fault dip angle for specified fault dip vector.
        * @param u array {u1,u2,u3} of components of fault dip vector.
        * @return fault dip angle, in degrees.
        */
        public static float faultDipFromDipVector(float[] u)
        {
            return faultDipFromDipVector(u[0], u[1], u[2]);
        }

        /**
        * Returns fault strike angle for specified fault strike vector.
        * @param v1 1st component of fault strike vector.
        * @param v2 2nd component of fault strike vector.
        * @param v3 3rd component of fault strike vector.
        * @return fault strike angle, in degrees.
        */
        public static float faultStrikeFromStrikeVector(
        float v1, float v2, float v3)
        {
            return range360(ArrayMath.toDegrees(Math.Atan2(v2, v3)));
        }

        /**
        * Returns fault strike angle for specified fault strike vector.
        * @param v array {v1,v2,v3} of components of fault strike vector.
        * @return fault strike angle, in degrees.
        */
        public static float faultStrikeFromStrikeVector(float[] v)
        {
            return faultStrikeFromStrikeVector(v[0], v[1], v[2]);
        }

        /**
        * Returns fault strike angle for specified fault normal vector.
        * The components w2 and w3 must not both be zero; that is, the
        * fault plane cannot be horizontal.
        * @param w1 1st component of fault normal vector.
        * @param w2 2nd component of fault normal vector.
        * @param w3 3rd component of fault normal vector.
        * @return fault strike angle, in degrees.
        */
        public static float faultStrikeFromNormalVector(
        float w1, float w2, float w3)
        {
            //Check.argument(w2 != 0.0f || w3 != 0.0f, "normal vector is not vertical");
            return range360(ArrayMath.toDegrees(Math.Atan2(-w3, w2)));
        }

        /**
        * Returns fault strike angle for specified fault normal vector.
        * The components w2 and w3 must not both be zero; that is, the
        * fault plane cannot be horizontal.
        * @param w array {w1,w2,w3} of components of fault normal vector.
        * @return fault strike angle, in degrees.
        */
        public static float faultStrikeFromNormalVector(float[] w)
        {
            return faultStrikeFromNormalVector(w[0], w[1], w[2]);
        }

        /**
        * Returns fault dip angle for specified fault normal vector.
        * @param w1 1st component of fault normal vector.
        * @param w2 2nd component of fault normal vector.
        * @param w3 3rd component of fault normal vector.
        * @return fault dip angle, in degrees.
        */
        public static float faultDipFromNormalVector(float w1, float w2, float w3)
        {
            return ArrayMath.toDegrees((float)Math.Acos(-w1));
        }

        /**
        * Returns fault dip angle for specified fault normal vector.
        * @param w array {w1,w2,w3} of components of fault normal vector.
        * @return fault dip angle, in degrees.
        */
        public static float faultDipFromNormalVector(float[] w)
        {
            return faultDipFromNormalVector(w[0], w[1], w[2]);
        }

        /**
        * Returns the cross-product of two specified vectors.
        * @param u the 1st vector.
        * @param v the 2nd vector.
        * @return array {w1,w2,w3} of components for vector w = u x v.
        */
        public static float[] crossProduct(float[] u, float[] v)
        {
            float u1 = u[0], u2 = u[1], u3 = u[2];
            float v1 = v[0], v2 = v[1], v3 = v[2];
            float w1 = u3 * v2 - u2 * v3;
            float w2 = u1 * v3 - u3 * v1;
            float w3 = u2 * v1 - u1 * v2;
            return new float[] { w1, w2, w3 };
        }

        /**
        * Returns angle in range [0,360] degrees.
        * @param phi angle, in degrees.
        * @return angle in range [0,360] degrees.
        */
        public static float range360(double phi)
        {
            while (phi < 0.0)
                phi += 360.0;
            while (phi >= 360.0)
                phi -= 360.0;
            return (float)phi;
        }

        /**
        * Returns angle in range [-180,180] degrees.
        * @param phi angle.
        * @return angle in range [-180,180] degrees.
        */
        public static float range180(double phi)
        {
            while (phi < -180.0)
                phi += 360.0;
            while (phi > 180.0)
                phi -= 360.0;
            return (float)phi;
        }

        ///////////////////////////////////////////////////////////////////////////
        // testing

        public static void main(String[] args)
        {
            float[] fp = {  0.0f, 90.0f,180.0f,270.0f,
        0.0f, 90.0f,180.0f,270.0f};
            float[] ft = { 90.0f, 90.0f, 90.0f, 90.0f,
        89.0f, 89.0f, 89.0f, 89.0f};
            for (int i = 0; i < fp.Length; ++i)
                test(fp[i], ft[i]);
        }
        public static void test(float phia, float thetaa)
        {
            float[] ua = faultNormalVectorFromStrikeAndDip(phia, thetaa);
            float phib = faultStrikeFromNormalVector(ua);
            float thetab = faultDipFromNormalVector(ua);
            float[] ub = faultNormalVectorFromStrikeAndDip(phib, thetab);
            assertEqual(ua, ub);
        }
        public static void assertEqual(float x, float y)
        {
            Debug.Assert(Math.Abs(x - y) < 0.01f);
        }
        public static void assertEqual(float[] x, float[] y)
        {
            for (int i = 0; i < x.Length; ++i)
                assertEqual(x[i], y[i]);
        }
        public static void trace(String s)
        {
            Console.WriteLine(s);
        }
    }
}
