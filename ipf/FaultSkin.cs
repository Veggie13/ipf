using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using edu.mines.jtk.awt;

namespace ipf
{
    public class FaultSkin : IEnumerable<FaultCell>
    {
        public IEnumerator<FaultCell> GetEnumerator()
        {
            foreach (FaultCell cell in _cellList)
            {
                yield return cell;
            }
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        private const long serialVersionUID = 1L;

        /**
        * Gets the cell that was the seed used to grow this skin.
        * @return the seed cell.
        */
        public FaultCell getSeed()
        {
            return _seed;
        }

        /**
        * Returns the number of cells in this skin.
        */
        public int size()
        {
            return _cellList.Count;
        }

        /**
        * Gets an array of cells in this skin.
        * @return array of cells.
        */
        public FaultCell[] getCells()
        {
            return _cellList.ToArray();
        }

        /**
        * Gets all cells in the specified skins.
        * @param skins array of skins for which to get cells.
        * @return array of cells.
        */
        public static FaultCell[] getCells(FaultSkin[] skins)
        {
            int ncell = countCells(skins);
            FaultCell[] cells = new FaultCell[ncell];
            int icell = 0;
            foreach (FaultSkin skin in skins)
                foreach (FaultCell cell in skin)
                    cells[icell++] = cell;
            return cells;
        }

        /**
        * Returns the total number of cells in the specified skins.
        * @param skins array of skins for which to count cells.
        * @return the total number of cells.
        */
        public static int countCells(FaultSkin[] skins)
        {
            int ncell = 0;
            foreach (FaultSkin skin in skins)
                ncell += skin.size();
            return ncell;
        }

        /**
        * Returns array of arrays of cells linked above and below.
        * @return array of arrays of linked cells; by reference, not by copy.
        */
        public FaultCell[][] getCellsAB()
        {
            if (_cellsAB != null)
                return _cellsAB;

            HashSet<FaultCell> cellSet = new HashSet<FaultCell>();
            List<FaultCell[]> cellsList = new List<FaultCell[]>();

            // For all cells in this skin, ...
            foreach (FaultCell cell in _cellList)
            {

                // If the cell is not already in an array, ...
                if (!cellSet.Contains(cell))
                {

                    // Search above for the top cell.
                    FaultCell c = cell;
                    for (FaultCell ca = c.ca; ca != null; ca = c.ca)
                        c = ca;

                    // Add the top cell and all cells below it to a new list.
                    List<FaultCell> cList = new List<FaultCell>();
                    for (; c != null; c = c.cb)
                    {
                        cList.Add(c);
                        cellSet.Add(c);
                    }

                    // Convert the list to an array and add it to the list of arrays.
                    cellsList.Add(cList.ToArray());
                }
            }
            Debug.Assert(_cellList.Count == cellSet.Count);

            // Convert the list of arrays to the array of arrays to be returned.
            _cellsAB = cellsList.ToArray();
            return _cellsAB;
        }

        /**
        * Returns array of arrays of cells linked left and right.
        * @return array of arrays of linked cells; by reference, not by copy.
        */
        public FaultCell[][] getCellsLR()
        {
            if (_cellsLR != null)
                return _cellsLR;

            HashSet<FaultCell> cellSet = new HashSet<FaultCell>();
            List<FaultCell[]> cellsList = new List<FaultCell[]>();

            // For all cells in this skin, ...
            foreach (FaultCell cell in _cellList)
            {

                // If the cell is not already in an array, ...
                if (!cellSet.Contains(cell))
                {

                    // Search left until we find no left nabor or until that left
                    // nabor is the cell with which we began the search.
                    FaultCell c = cell;
                    for (FaultCell cl = c.cl; cl != null && cl != cell; cl = c.cl)
                        c = cl;

                    // Remember the leftmost cell found and add it to a new list.
                    FaultCell cLeft = c;
                    List<FaultCell> cList = new List<FaultCell>();
                    cList.Add(c);
                    cellSet.Add(c);

                    // Add cells found to the right. Again beware of cycles.
                    for (c = c.cr; c != null && c != cLeft; c = c.cr)
                    {
                        cList.Add(c);
                        cellSet.Add(c);
                    }

                    // Convert the list to an array and add it to the list of arrays.
                    cellsList.Add(cList.ToArray());
                }
            }
            Debug.Assert(_cellList.Count == cellSet.Count);

            // Convert the list of arrays to the array of arrays to be returned.
            _cellsLR = cellsList.ToArray();
            checkCellArrays(_cellsLR);
            return _cellsLR;
        }

        /**
        * Smooths the normal vectors of cells in this skin.
        * @param nsmooth the number of smoothings.
        */
        public void smoothCellNormals(int nsmooth)
        {
            FaultCell.GetN getter = (FaultCell cell) => new float[] { cell.w1, cell.w2, cell.w3 };
            FaultCell.SetN setter = (FaultCell cell, float[] w) =>
            {
                float w1 = w[0];
                float w2 = w[1];
                float w3 = w[2];
                float ws = 1.0f / (float)Math.Sqrt(w1 * w1 + w2 * w2 + w3 * w3);
                w1 *= ws;
                w2 *= ws;
                w3 *= ws;
                cell.setNormalVector(w1, w2, w3);
            };
            for (int ismooth = 0; ismooth < nsmooth; ++ismooth)
                smoothN(getter, setter);
        }

        /**
        * Gets a cell nearest the centroid of this skin.
        * In illustrations, this cell is often a good representative.
        * @return the cell nearest the centroid.
        */
        public FaultCell getCellNearestCentroid()
        {
            float c1 = 0.0f;
            float c2 = 0.0f;
            float c3 = 0.0f;
            float cs = 0.0f;
            foreach (FaultCell c in _cellList)
            {
                c1 += c.fl * c.x1;
                c2 += c.fl * c.x2;
                c3 += c.fl * c.x3;
                cs += c.fl;
            }
            c1 /= cs;
            c2 /= cs;
            c3 /= cs;
            float dmin = float.MaxValue;
            FaultCell cmin = null;
            foreach (FaultCell c in _cellList)
            {
                float d = c.distanceSquaredTo(c1, c2, c3);
                if (d < dmin)
                {
                    cmin = c;
                    dmin = d;
                }
            }
            return cmin;
        }

        /**
        * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
        * these arrays, cells in this skin are represented by quads with specified
        * size, and colors corresponding to fault likelihoods.
        * @param size size (in samples) of the quads.
        * @param cmap colormap used to compute rgb colors from cell properties.
        * @param lhc true, if left-handed coordinate system; false, otherwise.
        */
        public float[][] getCellXyzUvwRgbForLikelihood(
        float size, ColorMap cmap, bool lhc)
        {
            return FaultCell.getXyzUvwRgbForLikelihood(size, cmap, getCells(), lhc);
        }

        /**
        * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
        * these arrays, cells in this skin are represented by quads with specified
        * size, and colors corresponding to fault throws.
        * @param size size (in samples) of the quads.
        * @param cmap colormap used to compute rgb colors from cell properties.
        * @param lhc true, if left-handed coordinate system; false, otherwise.
        */
        public float[][] getCellXyzUvwRgbForThrow(
        float size, ColorMap cmap, bool lhc)
        {
            return FaultCell.getXyzUvwRgbForThrow(size, cmap, getCells(), lhc);
        }

        /**
        * Gets arrays of packed cell coordinates for cell links.
        * Each returned array contains packed (x,y,z) cell coordinates for
        * exactly one above-below or left-right linked list of cells.
        * @return array of arrays of packed xyz cell coordinates.
        */
        public float[][] getCellLinksXyz()
        {
            FaultCell[][] cellsAB = getCellsAB();
            FaultCell[][] cellsLR = getCellsLR();
            int nsAB = cellsAB.Length; // number of segments for AB links
            int nsLR = cellsLR.Length; // number of segments for LR links
            int ns = nsAB + nsLR; // total number of segments
            float[][] xyz = new float[ns][];
            float[][] rgb = new float[ns][];
            for (int IS = 0; IS < ns; ++IS)
            { // for all segments, ...
                FaultCell[] cells = (IS < nsAB) ? cellsAB[IS] : cellsLR[IS - nsAB]; // the cells
                int ncell = cells.Length; // number of cells in thIS segment
                int np = ncell; // number of points in this segment
                if (IS >= nsAB && cells[0].cl == cells[ncell - 1]) // if a LR cycle, ...
                    ++np; // then add one more so we end with the starting point
                float[] xyzi = new float[3 * np]; // xyz for this segment
                for (int ip = 0, ic = 0; ip < np; ++ip)
                {
                    FaultCell cell = cells[ip % ncell]; // % to handle any LR cycle
                    xyzi[ic++] = cell.x3;
                    xyzi[ic++] = cell.x2;
                    xyzi[ic++] = cell.x1;
                }
                xyz[IS] = xyzi;
            }
            return xyz;
        }

        /**
        * Returns a fault skin read from a file with specified name.
        * @param fileName the fault skin file name.
        * @return the fault skin.
        */
        public static FaultSkin readFromFile(String fileName)
        {
            throw new NotImplementedException();
#if false
        FaultSkin skin = new FaultSkin();
        try {
        ArrayInputStream ais = new ArrayInputStream(fileName);
        int ncell = ais.readInt();
        int i1seed = ais.readInt();
        int i2seed = ais.readInt();
        int i3seed = ais.readInt();
        List<FaultCell> cellList = new List<FaultCell>(ncell);
        for (int icell=0; icell<ncell; ++icell) {
        float x1 = ais.readFloat();
        float x2 = ais.readFloat();
        float x3 = ais.readFloat();
        float fl = ais.readFloat();
        float fp = ais.readFloat();
        float ft = ais.readFloat();
        FaultCell cell = new FaultCell(x1,x2,x3,fl,fp,ft);
        cell.skin = skin;
        cellList.add(cell);
        cell.s1 = ais.readFloat();
        cell.s2 = ais.readFloat();
        cell.s3 = ais.readFloat();
        }
        FaultCell[] cells = cellList.toArray(new FaultCell[0]);
        FaultCellGrid fcg = new FaultCellGrid(cells);
        for (FaultCell cell:cellList) {
        int i1a = ais.readInt();
        int i2a = ais.readInt();
        int i3a = ais.readInt();
        int i1b = ais.readInt();
        int i2b = ais.readInt();
        int i3b = ais.readInt();
        int i1l = ais.readInt();
        int i2l = ais.readInt();
        int i3l = ais.readInt();
        int i1r = ais.readInt();
        int i2r = ais.readInt();
        int i3r = ais.readInt();
        if (i1a!=INULL) 
        cell.ca = fcg.get(i1a,i2a,i3a);
        if (i1b!=INULL) 
        cell.cb = fcg.get(i1b,i2b,i3b);
        if (i1l!=INULL) 
        cell.cl = fcg.get(i1l,i2l,i3l);
        if (i1r!=INULL) 
        cell.cr = fcg.get(i1r,i2r,i3r);
        }
        ais.close();
        skin._seed = fcg.get(i1seed,i2seed,i3seed);
        skin._cellList = cellList;
        } catch (Exception e) {
        throw new RuntimeException(e);
        }
        return skin;
#endif
        }

        /**
        * Writes a fault skin to a file with specified name.
        * @param fileName the fault skin file name.
        * @param skin the fault skin.
        */
        public static void writeToFile(String fileName, FaultSkin skin)
        {
            throw new NotImplementedException();
#if false
        try {
        ArrayOutputStream aos = new ArrayOutputStream(fileName);
        aos.writeInt(skin.size());
        aos.writeInt(skin._seed.i1);
        aos.writeInt(skin._seed.i2);
        aos.writeInt(skin._seed.i3);
        for (FaultCell cell:skin) {
        aos.writeFloat(cell.x1); 
        aos.writeFloat(cell.x2); 
        aos.writeFloat(cell.x3);
        aos.writeFloat(cell.fl); 
        aos.writeFloat(cell.fp); 
        aos.writeFloat(cell.ft);
        aos.writeFloat(cell.s1); 
        aos.writeFloat(cell.s2); 
        aos.writeFloat(cell.s3);
        }
        for (FaultCell cell:skin) {
        FaultCell[] nabors = new FaultCell[]{cell.ca,cell.cb,cell.cl,cell.cr};
        for (FaultCell nabor:nabors) {
        if (nabor!=null) {
        aos.writeInt(nabor.i1);
        aos.writeInt(nabor.i2);
        aos.writeInt(nabor.i3);
        } else {
        aos.writeInt(INULL);
        aos.writeInt(INULL);
        aos.writeInt(INULL);
        }
        }
        }
        aos.close();
        } catch (Exception e) {
        throw new RuntimeException(e);
        }
#endif
        }

        public static FaultSkin readFromFileSlow(String fileName)
        {
            throw new NotImplementedException();
#if false
            try
            {
                FileInputStream fis = new FileInputStream(fileName);
                ObjectInputStream ois = new ObjectInputStream(fis);
                FaultSkin skin = (FaultSkin)ois.readObject();
                ois.close();
                return skin;
            }
            catch (Exception e)
            {
                throw new RuntimeException(e);
            }
#endif
        }

        public static void writeToFileSlow(String fileName, FaultSkin skin)
        {
            throw new NotImplementedException();
#if false
            try
            {
                FileOutputStream fos = new FileOutputStream(fileName);
                ObjectOutputStream oos = new ObjectOutputStream(fos);
                oos.writeObject(skin);
                oos.close();
            }
            catch (Exception e)
            {
                throw new RuntimeException(e);
            }
#endif
        }

        /////////////////////////////////////////////////////////////////////////
        // package

        /**
        * Constructs an empty skin.
        */
        public FaultSkin()
        {
            _cellList = new List<FaultCell>();
        }

        /**
        * Adds the specified skinless cell to this skin.
        * @param cell the cell to be added.
        */
        void add(FaultCell cell)
        {
            Debug.Assert(cell.skin == null);
            cell.skin = this;
            if (_seed == null)
                _seed = cell;
            _cellList.Add(cell);
            _cellsAB = null;
            _cellsLR = null;
        }

        /**
        * Smooths one value stored in the cells of this skin. The value smoothed is
        * that accessed by the specified getter and setter. Each smoothed value is
        * an average of the values in a cell and its cell nabors. 
        */
        void smooth1(FaultCell.Get1 getter, FaultCell.Set1 setter)
        {
            int ncell = size();
            float[] vals = new float[ncell];
            float[] cnts = new float[ncell];
            FaultCell[] cellNabors = new FaultCell[4];
            for (int icell = 0; icell < ncell; ++icell)
            {
                FaultCell cell = _cellList[icell];
                float valCell = getter(cell);
                cellNabors[0] = cell.ca;
                cellNabors[1] = cell.cb;
                cellNabors[2] = cell.cl;
                cellNabors[3] = cell.cr;
                foreach (FaultCell cellNabor in cellNabors)
                {
                    if (cellNabor != null)
                    {
                        float valNabor = getter(cellNabor);
                        vals[icell] += valCell + valNabor;
                        cnts[icell] += 2.0f;
                    }
                }
            }
            for (int icell = 0; icell < ncell; ++icell)
            {
                FaultCell cell = _cellList[icell];
                float cnti = cnts[icell];
                float vali = vals[icell] / (cnti > 0.0f ? cnti : 1.0f);
                setter(cell, vali);
            }
        }

        /**
        * Smooths multiple values stored in the cells of this skin. The values
        * smoothed are those accessed by the specified getter and setter. Each
        * smoothed value is an average of the values in a cell and its cell nabors. 
        */
        void smoothN(FaultCell.GetN getter, FaultCell.SetN setter)
        {
            int ncell = size();
            int nval = getter(_seed).Length;
            float[][] vals = new float[ncell][];
            for (int xxx = 0; xxx < ncell; xxx++)
                vals[xxx] = new float[nval];
            float[] cnts = new float[ncell];
            FaultCell[] cellNabors = new FaultCell[4];
            for (int icell = 0; icell < ncell; ++icell)
            {
                FaultCell cell = _cellList[icell];
                float[] valsCell = getter(cell);
                cellNabors[0] = cell.ca;
                cellNabors[1] = cell.cb;
                cellNabors[2] = cell.cl;
                cellNabors[3] = cell.cr;
                foreach (FaultCell cellNabor in cellNabors)
                {
                    if (cellNabor != null)
                    {
                        float[] valsNabor = getter(cellNabor);
                        for (int ival = 0; ival < nval; ++ival)
                            vals[icell][ival] += valsCell[ival] + valsNabor[ival];
                        cnts[icell] += 2.0f;
                    }
                }
            }
            for (int icell = 0; icell < ncell; ++icell)
            {
                FaultCell cell = _cellList[icell];
                float cnti = cnts[icell];
                float scli = 1.0f / (cnti > 0.0f ? cnti : 1.0f);
                for (int ival = 0; ival < nval; ++ival)
                    vals[icell][ival] *= scli;
                setter(cell, vals[icell]);
            }
        }

        /////////////////////////////////////////////////////////////////////////
        // private

        private const int INULL = -int.MaxValue; // null index

        private FaultCell _seed; // cell in this skin with highest fl; null, if empty
        private List<FaultCell> _cellList; // list of cells in this skin
        private FaultCell[][] _cellsAB; // arrays of cells from above to below
        private FaultCell[][] _cellsLR; // arrays of cells from left to right

        private void checkCellArrays()
        {
            if (_cellsAB != null)
                checkCellArrays(_cellsAB);
            if (_cellsLR != null)
                checkCellArrays(_cellsLR);
        }

        private static void checkCellArrays(FaultCell[][] cells)
        {
            HashSet<FaultCell> cellSet = new HashSet<FaultCell>();
            int ncell = cells.Length;
            for (int icell = 0; icell < ncell; ++icell)
            {
                int mcell = cells[icell].Length;
                for (int jcell = 0; jcell < mcell; ++jcell)
                {
                    FaultCell c = cells[icell][jcell];
                    Debug.Assert(cellSet.Add(c));
                }
            }
        }

        private static void trace(String s)
        {
            Console.WriteLine(s);
        }
    }
}
