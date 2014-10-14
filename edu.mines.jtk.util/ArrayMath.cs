using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace edu.mines.jtk.util
{
    public static class ArrayMath
    {
        public static double toRadians(double angdeg)
        {
            return angdeg * Math.PI / 180;
        }

        public static double toDegrees(double angrad)
        {
            return angrad * 180 / Math.PI;
        }

        public static float toRadians(float angdeg)
        {
            return (float)toRadians((double)angdeg);
        }

        public static float toDegrees(float angrad)
        {
            return (float)toDegrees((double)angrad);
        }
    }
}
