﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using edu.mines.jtk.util;

namespace ipf
{
    public class ClosestPointTransform
    {
        /**
        * Constructs a closest-point transform.
        */
        public ClosestPointTransform()
        {
            _s1 = (float)1.0;
            _s2 = (float)1.0;
            _s3 = (float)1.0;
        }

        /**
        * Constructs a closest-point transform with specified scale factor.
        * @param s scale factor for all dimensions.
        */
        public ClosestPointTransform(double s)
        {
            _s1 = (float)s;
            _s2 = (float)s;
            _s3 = (float)s;
        }

        /**
        * Constructs a closest-point transform with specified scale factors.
        * @param s1 scale factor for 1st dimension.
        * @param s2 scale factor for 2nd and higher dimensions.
        */
        public ClosestPointTransform(double s1, double s2)
        {
            _s1 = (float)s1;
            _s2 = (float)s2;
            _s3 = (float)s2;
        }

        /**
        * Constructs a closest-point transform with specified scale factors.
        * @param s1 scale factor for 1st dimension.
        * @param s2 scale factor for 2nd dimension.
        * @param s3 scale factor for 3rd dimension.
        */
        public ClosestPointTransform(double s1, double s2, double s3)
        {
            _s1 = (float)s1;
            _s2 = (float)s2;
            _s3 = (float)s3;
        }

        /**
        * Applies this transform for a specified 1D sequence.
        * @param fnull value used to mark samples in f that are unknown.
        * @param f input array of values with at least one value known.
        * @param d output of distances to the nearest known value.
        * @param i output array of indices of nearest known value.
        */
        public void apply(float fnull, float[] f, float[] d, short[] i)
        {
            int n = f.Length;
            float[] g = new float[n];
            for (int j = 0; j < n; ++j)
                g[j] = (f[j] == fnull) ? HUGE : 0.0f;
            dt(_s1, g, d, i);
            for (int j = 0; j < n; ++j)
                d[j] = (float)Math.Sqrt(d[j]);
        }

        /**
        * Applies this transform for a specified 2D image.
        * @param fnull value used to mark samples in f that are unknown.
        * @param f input array of values with at least one value known.
        * @param d output of distances to the nearest known value.
        * @param i1 output array of 1st indices of nearest known value.
        * @param i2 output array of 2nd indices of nearest known value.
        */
        public void apply(
        float fnull, float[][] f, float[][] d,
        short[][] i1, short[][] i2)
        {
            int n1 = f[0].Length;
            int n2 = f.Length;
            Parallel.loop(n2, (int j2) =>
            {
                float[] f1 = new float[n1];
                float[] d1 = new float[n1];
                short[] k1 = new short[n1];
                for (int j1 = 0; j1 < n1; ++j1)
                    f1[j1] = (f[j2][j1] == fnull) ? HUGE : 0.0f;
                dt(_s1, f1, d1, k1);
                for (int j1 = 0; j1 < n1; ++j1)
                {
                    d[j2][j1] = d1[j1];
                    i1[j2][j1] = k1[j1];
                }
            });
            Parallel.loop(n1, (int j1) =>
            {
                int[] ij = new int[n2];
                float[] f2 = new float[n2];
                float[] d2 = new float[n2];
                short[] k1 = new short[n2];
                short[] k2 = new short[n2];
                for (int j2 = 0; j2 < n2; ++j2)
                    f2[j2] = d[j2][j1];
                dt(_s2, f2, d2, k2);
                for (int j2 = 0; j2 < n2; ++j2)
                {
                    d[j2][j1] = (float)Math.Sqrt(d2[j2]);
                    k1[j2] = i1[k2[j2]][j1];
                }
                for (int j2 = 0; j2 < n2; ++j2)
                {
                    i1[j2][j1] = k1[j2];
                    i2[j2][j1] = k2[j2];
                }
            });
        }

        /**
        * Applies this transform for a specified 3D image.
        * @param fnull value used to mark samples in f that are unknown.
        * @param f input array of values with at least one value known.
        * @param d output of distances to the nearest known value.
        * @param i1 output array of 1st indices of nearest known value.
        * @param i2 output array of 2nd indices of nearest known value.
        * @param i3 output array of 3rd indices of nearest known value.
        */
        public void apply(
        float fnull, float[][][] f, float[][][] d,
        short[][][] i1, short[][][] i2, short[][][] i3)
        {
            int n1 = f[0][0].Length;
            int n2 = f[0].Length;
            int n3 = f.Length;
            Parallel.loop(n3, (int j3) =>
            {
                { // axis 1
                    float[] f1 = new float[n1];
                    float[] d1 = new float[n1];
                    short[] k1 = new short[n1];
                    for (int j2 = 0; j2 < n2; ++j2)
                    {
                        for (int j1 = 0; j1 < n1; ++j1)
                            f1[j1] = (f[j3][j2][j1] == fnull) ? HUGE : 0.0f;
                        dt(_s1, f1, d1, k1);
                        for (int j1 = 0; j1 < n1; ++j1)
                        {
                            d[j3][j2][j1] = d1[j1];
                            i1[j3][j2][j1] = k1[j1];
                        }
                    }
                }
                { // axis 2
                    float[] f2 = new float[n2];
                    float[] d2 = new float[n2];
                    short[] k1 = new short[n2];
                    short[] k2 = new short[n2];
                    for (int j1 = 0; j1 < n1; ++j1)
                    {
                        for (int j2 = 0; j2 < n2; ++j2)
                            f2[j2] = d[j3][j2][j1];
                        dt(_s2, f2, d2, k2);
                        for (int j2 = 0; j2 < n2; ++j2)
                            k1[j2] = i1[j3][k2[j2]][j1];
                        for (int j2 = 0; j2 < n2; ++j2)
                        {
                            d[j3][j2][j1] = d2[j2];
                            i1[j3][j2][j1] = k1[j2];
                            i2[j3][j2][j1] = k2[j2];
                        }
                    }
                }
            });
            Parallel.loop(n2, (int j2) =>
            {
                // axis 3
                float[] f3 = new float[n3];
                float[] d3 = new float[n3];
                short[] k1 = new short[n3];
                short[] k2 = new short[n3];
                short[] k3 = new short[n3];
                for (int j1 = 0; j1 < n1; ++j1)
                {
                    for (int j3 = 0; j3 < n3; ++j3)
                        f3[j3] = d[j3][j2][j1];
                    dt(_s3, f3, d3, k3);
                    for (int j3 = 0; j3 < n3; ++j3)
                    {
                        k1[j3] = i1[k3[j3]][j2][j1];
                        k2[j3] = i2[k3[j3]][j2][j1];
                    }
                    for (int j3 = 0; j3 < n3; ++j3)
                    {
                        d[j3][j2][j1] = (float)Math.Sqrt(d3[j3]);
                        i1[j3][j2][j1] = k1[j3];
                        i2[j3][j2][j1] = k2[j3];
                        i3[j3][j2][j1] = k3[j3];
                    }
                }
            });
        }

        ///////////////////////////////////////////////////////////////////////////
        // private

        private float _s1, _s2, _s3;

        private static float HUGE = 1.0e20f;

        /**
        * Algorithm DT(f) of Felzenszwalb and Huttenlocher (2004),
        * augmented to compute the index i of the nearest sample.
        * @param s the distance between two adjacent samples.
        * @param f input array of accumulated distances squared.
        * @param d output array of distances squared to nearest samples.
        * @param i output array of indices of nearest samples.
        */
        private static void dt(float s, float[] f, float[] d, short[] i)
        {
            int n = f.Length;
            int[] v = new int[n];
            float[] z = new float[n + 1];
            z[0] = -HUGE;
            z[1] = HUGE;
            for (int q = 1, k = 0; q < n; ++q)
            {
                float r = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
                while (r <= z[k])
                {
                    --k;
                    r = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
                }
                ++k;
                v[k] = q;
                z[k] = r;
                z[k + 1] = HUGE;
            }
            for (int q = 0, k = 0; q < n; ++q)
            {
                while (z[k + 1] < q)
                    ++k;
                float qv = s * (q - v[k]);
                d[q] = qv * qv + f[v[k]];
                i[q] = (short)v[k];
            }
        }
    }
}
