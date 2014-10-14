using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace edu.mines.jtk.bench
{
    public class ArrayListBench
    {
        private const int INITIAL_CAPACITY = 8;

        public class FloatList
        {
            public int n;
            public float[] a = new float[INITIAL_CAPACITY];
            public void add(float f)
            {
                if (n == a.Length)
                {
                    float[] t = new float[2 * a.Length];
                    Array.Copy(a, 0, t, 0, n);
                    a = t;
                }
                a[n++] = f;
            }
            public float[] trim()
            {
                float[] t = new float[n];
                Array.Copy(a, 0, t, 0, n);
                return t;
            }
        }
    }
}
