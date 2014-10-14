using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace edu.mines.jtk.util
{
    public static class Parallel
    {
        public delegate void LoopInt(int i);

        public static void loop(int end, LoopInt body)
        {
            loop(0, end, 1, 1, body);
        }

        public static void loop(int begin, int end, int step, int chunk, LoopInt body)
        {
            if (!checkArgs(begin, end, step, chunk))
            {
                throw new ArgumentException();
            }

            var iS = Enumerable.Range(0, (end - begin) / step)
                .Select(i => begin + step * i);

            System.Threading.Tasks.Parallel.ForEach(iS, new Action<int>(body));
        }

        private static bool checkArgs(int begin, int end, int step, int chunk)
        {
            if (begin >= end)
            {
                return false;
            }
            if (step <= 0)
            {
                return false;
            }
            if (chunk <= 0)
            {
                return false;
            }

            return true;
        }
    }
}
