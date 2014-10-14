using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace edu.mines.jtk.dsp
{
    public class Sampling
    {
        /**
        * A default fraction used to test for equivalent sample values. 
        * By default, if the difference between two sample values does not 
        * exceed this fraction of the sampling interval, then those values 
        * are considered equivalent. This default is used when a tolerance 
        * is not specified explicitly when a sampling is constructed.
        */
        public const double DEFAULT_TOLERANCE = 1.0e-6;

        /**
        * Constructs a uniform sampling with specified count. The sampling 
        * interval (delta) is 1.0, and the first sample value is 0.0.
        * @param n the number (count) of samples; must be positive.
        */
        public Sampling(int n)
        {
            init(n, 1.0, 0.0, DEFAULT_TOLERANCE);
        }

        private void init(int n, double d, double f, double t)
        {
            //Check.argument(n>0,"n>0");
            //Check.argument(d>0.0,"d>0.0");
            _n = n;
            _d = d;
            _f = f;
            _v = null;
            _t = t;
            _td = _t * _d;
        }

        /**
        * Constructs a uniform sampling with specified parameters.
        * @param n the number (count) of samples; must be positive.
        * @param d the sampling interval (delta); must be positive.
        * @param f the first sample value.
        */
        public Sampling(int n, double d, double f)
        {
            init(n, d, f, DEFAULT_TOLERANCE);
        }

        /**
        * Constructs a sampling with specified parameters.
        * @param n the number (count) of samples; must be positive.
        * @param d the sampling interval (delta); must be positive.
        * @param f the first sample value.
        * @param t the sampling tolerance, expressed as fraction of delta.
        */
        public Sampling(int n, double d, double f, double t)
        {
            init(n, d, f, t);
        }

        /**
        * Constructs a sampling from the specified array of values. The values 
        * must be strictly increasing.
        * <p>
        * The constructed sampling may or may not be uniform, depending on the 
        * specified values and default sampling tolerance. If uniform (to within 
        * the default tolerance), then the array of values is discarded, and the 
        * sampling is represented by the count, sampling interval, and first 
        * sample value.
        * @param v the array of sampling values; must have non-zero length.
        */
        public Sampling(double[] v)
        {
            init(v, DEFAULT_TOLERANCE);
        }

        /**
        * Constructs a sampling from the specified array of values and tolerance.
        * The values must be strictly increasing. 
        * <p>
        * The constructed sampling may or may not be uniform, depending on the 
        * specified values and tolerance. If uniform (to within the specified
        * tolerance), then the array of values is discarded, and the sampling is 
        * represented by the count, sampling interval, and first sample value.
        * @param v the array of sampling values; must have non-zero length.
        * @param t the sampling tolerance, expressed as fraction of delta.
        */
        public Sampling(double[] v, double t)
        {
            init(v, t);
        }

        private void init(double[] v, double t)
        {
            //Check.argument(v.length>0,"v.length>0");
            //Check.argument(isIncreasing(v),"v is increasing");
            _n = v.Length;
            _d = (_n < 2) ? 1.0 : (v[_n - 1] - v[0]) / (_n - 1);
            _f = v[0];
            _t = t;
            _td = _t * _d;
            bool uniform = true;
            for (int i = 0; i < _n && uniform; ++i)
            {
                double vi = _f + i * _d;
                if (!almostEqual(v[i], vi, _td))
                    uniform = false;
            }
            _v = (uniform) ? null : v.ToList();
        }

        /**
        * Gets the number of samples.
        * @return the number of samples.
        */
        public int getCount()
        {
            return _n;
        }

        /**
        * Gets the sampling interval. If not uniformly sampled, the sampling
        * interval is the average difference between sample values.
        * @return the sampling interval; 1.0, if fewer than two samples.
        */
        public double getDelta()
        {
            return _d;
        }

        /**
        * Gets the first sample value.
        * @return the first sample value.
        */
        public double getFirst()
        {
            return _f;
        }

        /**
        * Gets the last sample value.
        * @return the last sample value.
        */
        public double getLast()
        {
            return (_v != null) ? _v[_n - 1] : _f + (_n - 1) * _d;
        }

        /**
        * Gets the sample value with specified index.
        * @param i the index.
        * @return the sample value.
        */
        public double getValue(int i)
        {
            //Check.index(_n,i);
            return value(i);
        }

        /**
        * Gets the sample values. If this sampling was constructed with an array 
        * of sample values, then the returned values are equivalent (equal to 
        * within the sampling tolerance) to the values in that array.
        * @return the sample values; returned by copy, not by reference.
        */
        public double[] getValues()
        {
            double[] v;
            if (_v != null)
            {
                v = _v.ToArray();
            }
            else
            {
                v = new double[_n];
                for (int i = 0; i < _n; ++i)
                    v[i] = _f + i * _d;
            }
            return v;
        }

        /**
        * Determines whether this sampling is uniform. A sampling is uniform 
        * if its values can be computed, to within the sampling tolerance, by
        * the expression v = f+i*d, for sampling indices i = 0, 1, ..., n-1.
        * Samplings with only one sample are considered to be uniform.
        * <p>
        * Note that, by this definition, samplings constructed with an array 
        * of sample values may or may not be uniform.
        * @return true, if uniform; false, otherwise.
        */
        public bool isUniform()
        {
            return _v == null;
        }

        /**
        * Determines whether this sampling is equivalent to the specified sampling.
        * Two samplings are equivalent if each of their sample values differs by 
        * no more than the sampling tolerance.
        * @param s the sampling to compare to this sampling.
        * @return true, if equivalent; false, otherwise.
        */
        public bool isEquivalentTo(Sampling s)
        {
            Sampling t = this;
            if (t.isUniform() != s.isUniform())
                return false;
            if (t.isUniform())
            {
                if (t.getCount() != s.getCount())
                    return false;
                double tiny = tinyWith(s);
                double tf = t.getFirst();
                double tl = t.getLast();
                double sf = s.getFirst();
                double sl = s.getLast();
                return almostEqual(tf, sf, tiny) && almostEqual(tl, sl, tiny);
            }
            else
            {
                double tiny = tinyWith(s);
                for (int i = 0; i < _n; ++i)
                {
                    if (!almostEqual(_v[i], s.value(i), tiny))
                        return false;
                }
                return true;
            }
        }

        /**
        * Determines whether this sampling is compatible with the specified sampling.
        * Two samplings are incompatible if their ranges of sample values overlap, 
        * but not all values in the overlapping parts are equivalent. Otherwise,
        * they are compatible.
        * @param s the sampling to compare to this sampling.
        * @return true, if compatible; false, otherwise.
        */
        public bool isCompatible(Sampling s)
        {
            Sampling t = this;
            int nt = t.getCount();
            int ns = s.getCount();
            double tf = t.getFirst();
            double sf = s.getFirst();
            double tl = t.getLast();
            double sl = s.getLast();
            int it = 0;
            int IS = 0;
            int jt = nt - 1;
            int js = ns - 1;
            if (tl < sf)
            {
                return true;
            }
            else if (sl < tf)
            {
                return true;
            }
            else
            {
                if (tf < sf)
                {
                    it = t.indexOf(sf);
                }
                else
                {
                    IS = s.indexOf(tf);
                }
                if (it < 0 || IS < 0)
                    return false;
                if (tl < sl)
                {
                    js = s.indexOf(tl);
                }
                else
                {
                    jt = t.indexOf(sl);
                }
                if (jt < 0 || js < 0)
                    return false;
                int mt = 1 + jt - it;
                int ms = 1 + js - IS;
                if (mt != ms)
                    return false;
                if (!t.isUniform() || !s.isUniform())
                {
                    double tiny = tinyWith(s);
                    for (jt = it, js = IS; jt != mt; ++jt, ++js)
                    {
                        if (!almostEqual(t.value(jt), s.value(js), tiny))
                            return false;
                    }
                }
                return true;
            }
        }

        /**
        * Returns the index of the sample with specified value. If this 
        * sampling has a sample value that equals (to within the sampling 
        * tolerance) the specified value, then this method returns the 
        * index of that sample. Otherwise, this method returns -1.
        * @param x the value.
        * @return the index of the matching sample; -1, if none.
        */
        public int indexOf(double x)
        {
            int i = -1;
            if (isUniform())
            {
                int j = (int)Math.Round((x - _f) / _d);
                if (0 <= j && j < _n && almostEqual(x, _f + j * _d, _td))
                    i = j;
            }
            else
            {
                int j = _v.BinarySearch(x);
                if (0 <= j)
                {
                    i = j;
                }
                else
                {
                    j = -(j + 1);
                    if (j > 0 && almostEqual(x, _v[j - 1], _td))
                    {
                        i = j - 1;
                    }
                    else if (j < _n && almostEqual(x, _v[j], _td))
                    {
                        i = j;
                    }
                }
            }
            return i;
        }

        /**
        * Returns the index of the sample nearest to the specified value.
        * @param x the value.
        * @return the index of the nearest sample.
        */
        public int indexOfNearest(double x)
        {
            int i;
            if (isUniform())
            {
                i = (int)Math.Round((x - _f) / _d);
                if (i < 0)
                    i = 0;
                if (i >= _n)
                    i = _n - 1;
            }
            else
            {
                i = _v.BinarySearch(x);
                if (i < 0)
                {
                    i = -(i + 1);
                    if (i == _n)
                    {
                        i = _n - 1;
                    }
                    else if (i > 0 && Math.Abs(x - _v[i - 1]) < Math.Abs(x - _v[i]))
                    {
                        --i;
                    }
                }
            }
            return i;
        }

        /**
        * Returns the value of the sample nearest to the specified value.
        * @param x the value.
        * @return the value of the nearest sample.
        */
        public double valueOfNearest(double x)
        {
            return getValue(indexOfNearest(x));
        }

        /**
        * Determines whether the specified index is in the bounds of this sampling.
        * An index is in bounds if in the range [0,count-1] of the first and last
        * sample indices.
        * @param i the index.
        * @return true, if in bounds; false, otherwise.
        */
        public bool isInBounds(int i)
        {
            return 0 <= i && i < _n;
        }

        /**
        * Determines whether the specified value is in the bounds of this sampling.
        * A value is in bounds if in the range [first,last] defined by the first 
        * and last sample values.
        * @param x the value.
        * @return true, if in bounds; false, otherwise.
        */
        public bool isInBounds(double x)
        {
            return getFirst() <= x && x <= getLast();
        }

        /**
        * Determines whether the specified value is in the bounds of this sampling,
        * which is assumed to be uniform. A value is in bounds if in the range
        * [first-0.5*delta,last+0.5*delta] defined by the first and last sample 
        * values and the sampling interval delta. In effect, this method extends
        * the bounds of this sampling by one-half sample when testing the value.
        * @param x the value.
        * @return true, if in bounds; false, otherwise.
        */
        public bool isInBoundsExtended(double x)
        {
            //Check.state(isUniform(),"sampling is uniform");
            double dhalf = 0.5 * _d;
            return getFirst() - dhalf <= x && x <= getLast() + dhalf;
        }

        /**
        * Gets the value for the specified index, assuming uniform sampling.
        * The index and the returned value need not be in the bounds of this 
        * sampling, which must be uniform. That is, the specified index may be 
        * less than zero or greater than or equal to the number of samples.
        * @param i the index.
        * @return the value.
        */
        public double getValueExtended(int i)
        {
            //Check.state(isUniform(),"sampling is uniform");
            return _f + i * _d;
        }

        /**
        * Returns the index of the sample nearest the specified value, assuming 
        * uniform sampling. The value and the returned index need not be in the 
        * bounds of this sampling, which must be uniform. Specifically, the
        * returned index may be less than zero or greater than or equal to the 
        * number of samples.
        * @param x the value.
        */
        public int indexOfNearestExtended(double x)
        {
            //Check.state(isUniform(),"sampling is uniform");
            return (int)Math.Round((x - _f) / _d);
        }

        /**
        * Returns the value of the sample nearest to the specified value, 
        * assuming uniform sampling. The specified and returned values need 
        * not be in the bounds of this sampling, which must be uniform.
        * @param x the value.
        * @return the value of the nearest sample.
        */
        public double valueOfNearestExtended(double x)
        {
            return getValueExtended(indexOfNearestExtended(x));
        }

        /**
        * Returns the index of the floor of the specified value.
        * The floor is the largest sampled value less than or equal to the 
        * specified value. The value and the returned index need not be in 
        * the bounds of this sampling, which must be uniform.
        * Specifically, the returned index may be less than zero or greater 
        * than or equal to the number of samples.
        * @param x the value.
        */
        public int indexOfFloorExtended(double x)
        {
            //Check.state(isUniform(),"sampling is uniform");
            double xn = (x - _f) / _d;
            int ix = (int)xn;
            return (ix <= xn) ? ix : ix - 1;
        }

        /**
        * Returns the normalized difference between a specified value and 
        * the sampled value for a specified index. Normalized difference is
        * the difference (specified value minus sampled value) divided by 
        * the sampling interval. The specified value and index not be in 
        * the bounds of this sampling, which must be uniform.
        * @param x the value.
        * @param i the index.
        */
        public double normalizedDifferenceExtended(double x, int i)
        {
            //Check.state(isUniform(),"sampling is uniform");
            return (x - (_f + i * _d)) / _d;
        }

        /**
        * Determines the overlap between this sampling and the specified sampling.
        * Both the specified sampling and this sampling represent a first-to-last
        * range of sample values. The overlap is a contiguous set of samples that 
        * have values that are equal, to within the minimum sampling tolerance of 
        * the two samplings. This set is represented by an array of three ints, 
        * {n,it,is}, where n is the number of overlapping samples, and it and is 
        * denote the indices of the first samples in the overlapping parts of this 
        * sampling (t) and the specified sampling (s). There exist three cases.
        * <ul><li>
        * The ranges of sample values overlap, and all values in the overlapping 
        * parts are equivalent. In this case, the two samplings are compatible; 
        * and this method returns the array {n,it,is}, as described 
        * above.
        * </li><li>
        * The ranges of sample values in the two samplings do not overlap.
        * In this case, the two samplings are compatible; and this method 
        * returns either {0,nt,0} or {0,0,ns}, depending on whether all
        * sample values in this sampling are less than or greater than
        * those in the specified sampling, respectively.
        * </li><li>
        * The ranges of sample values overlap, but not all values in the 
        * overlapping parts are equivalent. In this case, the two samplings 
        * are incompatible; and this method returns null.
        * </li></ul>
        * @param s the sampling to compare with this sampling.
        * @return the array {n,it,is} that describes the overlap; null, if 
        *  the specified sampling is incompatible with this sampling.
        */
        public int[] overlapWith(Sampling s)
        {
            Sampling t = this;
            int nt = t.getCount();
            int ns = s.getCount();
            double tf = t.getFirst();
            double sf = s.getFirst();
            double tl = t.getLast();
            double sl = s.getLast();
            int it = 0;
            int IS = 0;
            int jt = nt - 1;
            int js = ns - 1;
            if (tl < sf)
            {
                return new int[] { 0, nt, 0 };
            }
            else if (sl < tf)
            {
                return new int[] { 0, 0, ns };
            }
            else
            {
                if (tf < sf)
                {
                    it = t.indexOf(sf);
                }
                else
                {
                    IS = s.indexOf(tf);
                }
                if (it < 0 || IS < 0)
                    return null;
                if (tl < sl)
                {
                    js = s.indexOf(tl);
                }
                else
                {
                    jt = t.indexOf(sl);
                }
                if (jt < 0 || js < 0)
                    return null;
                int mt = 1 + jt - it;
                int ms = 1 + js - IS;
                if (mt != ms)
                    return null;
                if (!t.isUniform() || !s.isUniform())
                {
                    double tiny = tinyWith(s);
                    for (jt = it, js = IS; jt != mt; ++jt, ++js)
                    {
                        if (!almostEqual(t.value(jt), s.value(js), tiny))
                            return null;
                    }
                }
                return new int[] { mt, it, IS };
            }
        }

        /**
        * Returns the union of this sampling with the specified sampling. This
        * union is possible if and only if the two samplings are compatible.
        * <p>
        * If the two samplings do not overlap, this method does not create
        * samples within any gap that may exist between them. In other words,
        * the number of samples in the sampling returned is exactly nt+ns-n, 
        * where nt is the number of samples in this sampling, ns is the number 
        * of samples in the specified sampling, and n is the number of samples 
        * with equivalent values in any overlapping parts of the two samplings.
        * If the samplings do not overlap, then n = 0. One consequence of this
        * behavior is that the union of two uniform samplings with the same 
        * sampling interval may be non-uniform.
        * <p>
        * This method returns a new sampling; it does not modify this sampling.
        * @see #overlapWith(Sampling)
        * @param s the sampling to merge with this sampling.
        * @return the merged sampling; null, if no merge is possible.
        */
        public Sampling mergeWith(Sampling s)
        {
            Sampling t = this;
            int[] overlap = t.overlapWith(s);
            if (overlap == null)
                return null;
            int n = overlap[0];
            int it = overlap[1];
            int IS = overlap[2];
            int nt = t.getCount();
            int ns = s.getCount();
            int nm = nt + ns - n;
            if (n > 0 && t.isUniform() && s.isUniform())
            {
                double dm = t.getDelta();
                double fm = (it == 0) ? s.getFirst() : t.getFirst();
                return new Sampling(nm, dm, fm);
            }
            else
            {
                double[] vm = new double[nm];
                int jm = 0;
                for (int jt = 0; jt < it; ++jt)
                    vm[jm++] = t.value(jt);
                for (int js = 0; js < IS; ++js)
                    vm[jm++] = s.value(js);
                for (int jt = it; jt < it + n; ++jt)
                    vm[jm++] = t.value(jt);
                for (int jt = it + n; jt < nt; ++jt)
                    vm[jm++] = t.value(jt);
                for (int js = IS + n; js < ns; ++js)
                    vm[jm++] = s.value(js);
                return new Sampling(vm);
            }
        }

        /**
        * Shifts this sampling.
        * <p>
        * This method returns a new sampling; it does not modify this sampling.
        * @param s the value (shift) to add to this sampling's values.
        * @return the new sampling.
        */
        public Sampling shift(double s)
        {
            if (_v == null)
            {
                return new Sampling(_n, _d, _f + s, _t);
            }
            else
            {
                double[] v = new double[_n];
                for (int i = 0; i < _n; ++i)
                    v[i] = _v[i] + s;
                return new Sampling(v, _t);
            }
        }

        /**
        * Prepends samples to this sampling.
        * If this sampling is not uniform, prepended sample values are computed 
        * to preserve the average difference between adjacent sample values.
        * <p>
        * This method returns a new sampling; it does not modify this sampling.
        * @param m the number of new samples prepended to this sampling.
        * @return the new sampling.
        */
        public Sampling prepend(int m)
        {
            int n = _n + m;
            double f = _f - m * _d;
            if (_v == null)
            {
                return new Sampling(n, _d, f, _t);
            }
            else
            {
                double[] v = new double[n];
                for (int i = 0; i < m; ++i)
                    v[i] = f + i * _d;
                for (int i = m; i < n; ++i)
                    v[i] = _v[i - m];
                return new Sampling(v, _t);
            }
        }

        /**
        * Appends samples to this sampling.
        * If this sampling is not uniform, appended sample values are computed 
        * to preserve the average difference between adjacent sample values.
        * <p>
        * This method returns a new sampling; it does not modify this sampling.
        * @param m the number of new samples appended to this sampling.
        * @return the new sampling.
        */
        public Sampling append(int m)
        {
            int n = _n + m;
            if (_v == null)
            {
                return new Sampling(n, _d, _f, _t);
            }
            else
            {
                double[] v = new double[n];
                for (int i = 0; i < _n; ++i)
                    v[i] = _v[i];
                for (int i = _n; i < n; ++i)
                    v[i] = _f + i * _d;
                return new Sampling(v, _t);
            }
        }

        /**
        * Decimates this sampling. Beginning with the first sample, keeps only 
        * every m'th sample, while discarding the others in this sampling. If 
        * this sampling has n values, the new sampling will have 1+(n-1)/m values.
        * <p>
        * This method returns a new sampling; it does not modify this sampling.
        * @param m the factor by which to decimate; must be positive.
        * @return the new sampling.
        */
        public Sampling decimate(int m)
        {
            //Check.argument(m>0,"m>0");
            int n = 1 + (_n - 1) / m;
            if (_v == null)
            {
                return new Sampling(n, m * _d, _f, _t);
            }
            else
            {
                double[] v = new double[n];
                for (int i = 0, j = 0; i < n; ++i, j += m)
                    v[i] = _v[j];
                return new Sampling(v, _t);
            }
        }

        /**
        * Interpolates this sampling. Inserts m-1 evenly spaced samples between 
        * each of the samples in this sampling. If this sampling has n values, 
        * the new sampling will have 1+(n-1)*m values.
        * <p>
        * This method returns a new sampling; it does not modify this sampling.
        * @param m the factor by which to interpolate.
        * @return the new sampling.
        */
        public Sampling interpolate(int m)
        {
            //Check.argument(m>0,"m>0");
            int n = _n + (_n - 1) * (m - 1);
            if (_v == null)
            {
                return new Sampling(n, _d / m, _f, _t);
            }
            else
            {
                double[] v = new double[n];
                v[0] = _v[0];
                for (int i = 1, j = m; i < _n; ++i, j += m)
                {
                    v[j] = _v[i];
                    double dk = (v[j] - v[j - m]) / m;
                    double vk = v[j - m];
                    for (int k = j - m + 1; k < j; ++k, vk += dk)
                        v[k] = vk;
                }
                return new Sampling(v, _t);
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        // private

        private int _n; // number of samples
        private double _d; // sampling interval
        private double _f; // value of first sample
        private List<double> _v; // array[n] of sample values; null, if uniform
        private double _t; // sampling tolerance, as a fraction of _d
        private double _td; // sampling tolerance _t multiplied by _d

        private double value(int i)
        {
            return (_v != null) ? _v[i] : _f + i * _d;
        }

        private bool almostEqual(double v1, double v2, double tiny)
        {
            double diff = v1 - v2;
            return (diff < 0.0) ? -diff < tiny : diff < tiny;
        }

        private double tinyWith(Sampling s)
        {
            return (_td < s._td) ? _td : s._td;
        }
    }
}
