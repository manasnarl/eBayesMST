using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using eBayesMSTClass;
using ILNumerics;
using System.Diagnostics;
using System.ComponentModel;

namespace eBayesMST

{
    /// <summary>
    /// Interaction logic for eBayesAlgorithm.xaml
    /// </summary>
    public partial class eBayesMSTAlgorithm : Window
    {
        eBayesMSTClass.eBayesMSTClass clseBayes = new eBayesMSTClass.eBayesMSTClass();
        public eBayesMSTAlgorithm()
        {
            InitializeComponent();
        }
        double[] R_sig = null;
       
        private Double[] MultiSignalWaveletDecomposition_eBayes(double[] R_sig, int nlvl, double[] LoD, double[] HiD, double[] LoR, double[] HiR)
        {
            //Algorithm Step 11.1 Decompose R_sig using Multi signal Wavelet Decomposition 
            int nfft = R_sig.Length;
            lstListofModules.Items.Add("Algorithm Step 11.1 Decompose R_sig using Multi signal Wavelet Decomposition completed ");

            //Algorithm Step 11.1.1. Extend, Decompose and Extract → xDec
            double[] WD = new double[nfft];
            double[] xDec_CA = new double[nfft];
            double[][] xDec_CD = new double[nfft][];
            Int32[] SHE = new Int32[nlvl + 2];
            SHE[SHE.Length - 1] = nfft;
            lstListofModules.Items.Add("Algorithm Step 11.1.1. Extend, Decompose and Extract → xDec completed");

            //Algorithm Step 11.1.1.1. Initialize Approximation Coefficients → xDec.Ca
            //Algorithm Step 11.1.1.2. Initialize Detailed Coefficients for nlvl levels → xDec.Cd
            mdwtdec(R_sig, nlvl, LoD, HiD, LoR, HiR, out xDec_CA, out xDec_CD, SHE);
            SHE[0] = xDec_CA.Length;

            //Algorithm Step 11.2. Calculate Normalization Factor → norf
            double norf = (double)ILMath.divide(1, -(ILMath.multiply(ILMath.sqrt((double)2), -0.476936276204470)));//-0.476936276204470=erfcinv(2*0.75)
            lstListofModules.Items.Add("Algorithm Step 11.2. Calculate Normalization Factor → norf completed");
            

            double[] vscale = null;

            //Algorithm Step 11.3. Calculate median of xDec.Cd at level1 → Med
            ILArray<double> Med = xDec_CD[0];
            lstListofModules.Items.Add("Algorithm Step 11.3. Calculate median of xDec.Cd at level1 → Med completed");

            //Algorithm Step 11.4. Store(norf * Med) → vScale
            vscale = ILMath.multiply(norf, ILMath.median(ILMath.abs(Med))).ToArray();
            lstListofModules.Items.Add("Algorithm Step 11.4. Store(norf * Med) → vScale completed");

            //Algorithm Step 11.5.Repeat calculation of eBayes-threshold for at all nlvl levels
            lstListofModules.Items.Add("Algorithm Step 11.5. Repeat calculation of eBayes-threshold for at all nlvl levels completed");
            for (int lev = 0; lev < nlvl; lev++)
            {
                //Algorithm Step 11.5.1 Using xDec.Cd at the current level and vScale → xDecL
                xDec_CD[lev] = ebayesThresh(xDec_CD[lev], vscale);
            }


            lstListofModules.Items.Add("Algorithm Step 11.5.1 Using xDec.Cd at the current level and vScale → xDecL completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.1. Using xDec.Cd at the current level and vScale → xDecL completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.2. Recalculate Normalization Factor for xDecL → xnorf completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.3. Calculate median for xDecL → mxDecL completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.4. Subtract mxDecL from xDecL to recalculate median → mxDecLm completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.5. Replicate mxDecLm * xnorf to the size of xDecL → rtStdEst completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.6. Divide xDecL with rtStdEst for unit standard deviation → x completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.7. Calculate prior weight from x → mu completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.8. Calculate Posterior median estimates until it converges or reaches maximum iterations using mu → muhat completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.9. Multiply muhat with rtStdEst → omuhat completed");
            lstListofModules.Items.Add("Algorithm Step 11.5.10. Replace xDec.Cd with current level omuhat completed");


            //Algorithm Step 11.6.3 F_o as Segment output
            double[] F_o = mdwtrec(LoR, HiR, xDec_CA, xDec_CD, nlvl, SHE);
            
           
            //return F_o as Segment output
            return F_o;

        }
        public double[] cauchymedzero(double[] muhat, double[] x, double[] weight)
        {
            double[] y = ILMath.subtract((ILArray<double>)x, (ILArray<double>)muhat).ToArray();
            double[] fx = clseBayes.gausspdf(y, 0, 1);
            double[] a = clseBayes.gausscdf(y, 0, 1, "lower");
            double[] minmuhat = ILMath.multiplyElem(muhat, (double)-1).ToArray();

            double[] b = clseBayes.gausscdf(minmuhat, 0, 1, "lower");
            double[] c = clseBayes.gausspdf(muhat, 0, 1);
            double[] yr = new double[x.Length];
            for (int i = 0; i < yr.Length; i++)
                yr[i] = a[i] - (x[i] * fx[i]) + ((((x[i] * muhat[i]) - 1) * fx[i] * b[i]) / c[i]);
            double[] temp = ((ILMath.multiplyElem(ILMath.pow((ILArray<double>)x, 2.0f), (ILMath.divide((double)1, weight) - 1.0f))) - 1.0f).ToArray();
            double[] y1 = ILMath.add(1, ILMath.multiplyElem(ILMath.exp(ILMath.divide(-ILMath.pow((ILArray<double>)x, 2), 2)), temp)).ToArray();
            double[] z = ILMath.subtract(ILMath.divide((ILArray<double>)y1, 2), yr).ToArray();

            return z;
        }
        public void intervalSolve(double[] zeromd, double lo, double hi, int maxiter, double[] magdata, double[] weight, out double[] muhat, out List<double> delta)
        {

            int len = zeromd.Length;
            int[] m = new int[len];
            double[] tlo = new double[len]; double[] thi = new double[len];
            for (int i = 0; i < len; i++)
            {
                tlo[i] = lo; thi[i] = hi;
            }
            double Tol = 1e-09;
            int numiter = 0;
            double conTol = double.PositiveInfinity;
            List<double> temp_delta = new List<double>();
            while (conTol > Tol)
            {
                numiter = numiter + 1;
                double[] ArrmidPoint = ILMath.divide(ILMath.add((ILArray<double>)tlo, thi), 2).ToArray();
                double[] fmidpoint = cauchymedzero(ArrmidPoint, magdata, weight);
                double[] idx = new double[len];
                for (int i = 0; i < len; i++)
                    if (fmidpoint[i] <= zeromd[i]) idx[i] = 1;
                    else idx[i] = 0;
                for (int i = 0; i < len; i++)
                    if (idx[i] == 1) tlo[i] = ArrmidPoint[i];
                    else thi[i] = ArrmidPoint[i];
                temp_delta.Add((double)ILMath.max(ILMath.abs(ILMath.subtract((ILArray<double>)thi, tlo))));
                double temp_max = temp_delta[numiter - 1];
                conTol = (temp_max);
                if (numiter > maxiter) break;
            }
            delta = temp_delta;
            muhat = ILMath.divide(ILMath.add((ILArray<double>)tlo, thi), 2).ToArray();
        }
        private void mdwtdec(double[] x, int level, ILArray<double> LoD, ILArray<double> HiD, ILArray<double> LoR, ILArray<double> HiR, out double[] xDec_CA, out double[][] xDec_CD, int[] SHE)
        {
            //Algorithm Step 11.1.2 Store the symmetric half-point extension → SHE
            
            
            int First = 2; xDec_CD = new double[level][]; xDec_CA = new double[x.Length];
            for (int j = level + 1; j > 1; j--)
            {
                int lf = LoD.Length;
                int lx = x.Length;
                int dcol = lf - 1;
                int lenEXT = lf - 1;
                int lenKEPT = lx + lf - 1;
                List<int> idxCOL = new List<int>();
                int Checking = 1;
                for (int i = 0; i < Checking; i++)
                {
                    if (i == 0)
                    {
                        idxCOL.Add(First + dcol); Checking++;
                    }
                    else
                    {
                        int d = idxCOL[i - 1] + 2;
                        if (d <= (lenKEPT + dcol))
                        {
                            idxCOL.Add(d);
                            Checking++;
                        }
                        else Checking = 0;
                    }

                }
                int[] I = new int[lenKEPT + dcol];
                bool zeroaccept = false; bool nfftAccept = false;
                for (int i = 0; i < I.Length; i++)
                {
                    if (i == 0) I[i] = lenEXT;
                    else
                    {
                        if (!zeroaccept)
                        {
                            if ((I[i - 1] - 1) > 0)
                                I[i] = I[i - 1] - 1;
                            else if ((I[i - 1] - 1) == 0)
                            {
                                I[i] = I[i - 1];
                                zeroaccept = true;
                            }
                        }
                        else if (!nfftAccept && zeroaccept)
                        {
                            if ((I[i - 1] + 1) < lx + 1)
                                I[i] = I[i - 1] + 1;
                            else if ((I[i - 1] + 1) == lx + 1)
                            {
                                I[i] = I[i - 1];
                                nfftAccept = true;
                            }
                        }
                        else
                            I[i] = I[i - 1] - 1;
                    }

                }

                double[] y = new double[I.Length];
                for (int i = 0; i < y.Length; i++)
                {
                    y[i] = x[I[i] - 1];
                }
                //Algorithm Step 11.1.3. Do the 2D convolution of SHE with LoD → xDec.Ca
                double[] conva = clseBayes.Convlution2D(y, LoD);
                double[] a = new double[idxCOL.Count];

                for (int i = 0; i < a.Length; i++)
                {
                    a[i] = conva[idxCOL[i] - 1];
                }
                

                //Algorithm Step 11.1.4. Do the 2D convolution of SHE with HiD → xDec.Cd
                double[] Convd = clseBayes.Convlution2D(y, HiD);
                double[] dd = new double[idxCOL.Count];

                for (int i = 0; i < dd.Length; i++)
                {
                    dd[i] = Convd[idxCOL[i] - 1];
                }

                x = a;
                xDec_CA = new double[x.Length];
                xDec_CA = a;
                xDec_CD[level + 1 - j] = new Double[dd.Length];
                for (int k = 0; k < dd.Length; k++)
                    xDec_CD[level + 1 - j][k] = dd[k];

                

                SHE[j - 1] = dd.Length;

            }
            lstListofModules.Items.Add("Algorithm Step 11.1.2 Store the symmetric half-point extension → SHE completed");
            lstListofModules.Items.Add("Algorithm Step 11.1.3. Do the 2D convolution of SHE with LoD → xDec.Ca completed");
            lstListofModules.Items.Add("Algorithm Step 11.1.4. Do the 2D convolution of SHE with HiD → xDec.Cd completed");
        }

        private double[] mdwtrec(ILArray<double> LoR, ILArray<double> HiR, double[] CA, double[][] CD, int level, int[] sx)
        {
            int levMin = 0;
            int levMax = level;
            double[] x = CA;
            Double[] F_o = null;
            //Algorithm Step 11.6 Multi signal Wavelet Reconstruction at all levels of xDec.Ca and Repeat the following for all nlvl levels from max-level to min - level
            
            for (int j = levMax; j > levMin; j--)
            {

                int p1 = levMax + 2 - j;
                int lenkept = sx[p1];
                int sx2 = 2 * x.Length;
                double[] y1 = new double[sx2 - 1];
                int k = 0;
                for (int i = 0; i < y1.Length; i++)
                {
                    if (i % 2 == 0) { y1[i] = x[k]; k++; }
                    else y1[i] = 0;
                }
                //Algorithm Step 11.6.1. Up - Sample, 2D convolution using xDec.Ca and LoR → lxDec
                double[] lxDec = clseBayes.Convlution2D(y1, LoR);
                int sy = lxDec.Length;
                if (lenkept > sy) lenkept = sy;
                double d = (sy - lenkept) / 2;
                int first = 1 + Convert.ToInt32(Math.Floor(d));
                int Last = sy - Convert.ToInt32(Math.Ceiling(d));

                double[] SumY1 = new double[(Last - first) + 1];

                int m = first - 1;
                for (int i = 0; i < SumY1.Length; i++)
                {
                    SumY1[i] = lxDec[m]; m++;
                }


                sx2 = 2 * CD[j - 1].Length;
                double[] y2 = new double[sx2 - 1];
                k = 0;
                for (int i = 0; i < y2.Length; i++)
                {
                    if (i % 2 == 0) { y2[i] = CD[j - 1][k]; k++; }
                    else y2[i] = 0;
                }
                //Algorithm Step 11.6.2. Up-Sample, 2D convolution using xDec.Ca and HiR → hxDec
                double[] hxDec = clseBayes.Convlution2D(y2, HiR);
                sy = hxDec.Length;
                if (lenkept > sy) lenkept = sy;
                d = (sy - lenkept) / 2;
                first = 1 + Convert.ToInt32(Math.Floor(d));
                Last = sy - Convert.ToInt32(Math.Ceiling(d));

                double[] SumY2 = new double[(Last - first) + 1];

                m = first - 1;
                for (int i = 0; i < SumY2.Length; i++)
                {
                    SumY2[i] = hxDec[m]; m++;
                }

                //Algorithm Step 11.6.3. Add lxDec and hxDec → F_o as Segment output
                F_o = new double[lenkept];

                for (int i = 0; i < F_o.Length; i++)
                    F_o[i] = SumY1[i] + SumY2[i];


                x = F_o;

            }
            lstListofModules.Items.Add("Algorithm Step 11.6 Multi signal Wavelet Reconstruction at all levels of xDec.Ca and Repeat the following for all nlvl levels from max-level to min - level completed");
            lstListofModules.Items.Add("Algorithm Step 11.6.1. Up - Sample, 2D convolution using xDec.Ca and LoR → lxDec completed");
            lstListofModules.Items.Add("Algorithm Step 11.6.2. Up-Sample, 2D convolution using xDec.Ca and HiR → hxDec completed");
            lstListofModules.Items.Add("Algorithm Step 11.6.3. Add lxDec and hxDec → F_o as Segment output completed");

            //return Segment output F_o
            return F_o;
        }

        private double[] ebayesThresh(double[] xDecL, double[] vscale)
        {
            //Algorithm Step 11.5.1. Using xDec.Cd at the current level and vScale → xDecL
            int maxiter = 50;
            double minstd = 1e-9;
            int m = xDecL.Length;
            double[] omuhat = xDecL;

            ILArray<double> temp = null;
            double[] rtStdEst = null;
            int[] size = new int[1] { m };

            //Algorithm Step 11.5.2. Recalculate Normalization Factor for xDecL → xnorf
            double xnorf = (double)ILMath.divide(1, -(ILMath.multiply(ILMath.sqrt((double)2), -0.476936276204470)));
            

            //Algorithm Step 11.5.3. Calculate median for xDecL → mxDecL
            ILArray<double> mxDecL = ILMath.median((ILArray<double>)xDecL);
            

            //Algorithm Step 11.5.4. Subtract mxDecL from xDecL to recalculate median → mxDecLm
            ILArray<double> mxDecLm = ILMath.median(ILMath.abs(ILMath.subtract((ILArray<double>)xDecL, mxDecL)));
            

            //Algorithm Step 11.5.5. Replicate mxDecLm * xnorf to the size of xDecL → rtStdEst
            temp = ILMath.multiplyElem(xnorf, mxDecLm).ToArray();
            rtStdEst = ILMath.repmat(temp, m).ToArray();

            for (int i = 0; i < rtStdEst.Length; i++)
                if (rtStdEst[i] < minstd)
                    rtStdEst[i] = minstd;

            

            //Algorithm Step 11.5.6. Divide xDecL with rtStdEst for unit standard deviation → x
            double[] x = ILMath.divide((ILArray<double>)xDecL, rtStdEst).ToArray();
            


            //Algorithm Step 11.5.7. Calculate prior weight from x → mu
            double mu = clseBayes.weightfromData(x, maxiter, "decimated");
            

            double[] muhat = new double[x.Length];

            //Algorithm Step 11.5.8. Calculate Posterior median estimates until it converges or reaches maximum iterations using mu → muhat
            muhat = postmedCauchy(x, mu, maxiter);
            

            //Algorithm Step 11.5.9. Multiply muhat with rtStdEst → omuhat
            omuhat = ILMath.multiplyElem((ILArray<double>)muhat, rtStdEst).ToArray();
            

            //Algorithm Step 11.5.10. Replace xDec.Cd with current level omuhat 
            return omuhat;
            
        }
        public double[] postmedCauchy(double[] data, double weight, int maxiter)
        {
            double[] muhat = new double[data.Length];
            int M = data.Length;
            int N = 1;
            int[] muhatlen = new int[M];

            double[] Weight = new double[M];
            for (int i = 0; i < M; i++)
                Weight[i] = weight;
            double[] magdata = ILMath.abs((ILArray<double>)data).ToArray();

            double[] magdatatmp = new double[M];
            for (int i = 0; i < M; i++)
                magdatatmp[i] = magdata[i];

            double[] idx = new double[M];
            for (int i = 0; i < M; i++)
                if (magdata[i] < 20) idx[i] = 1;
                else idx[i] = 0;
            for (int i = 0; i < M; i++)
                if (idx[i] == 0) magdata[i] = double.NaN;
            double lo = 0; List<double> delta = new List<double>();
            double[] zeromd = new double[magdata.Length];
            intervalSolve(zeromd, lo, magdata.Max(), maxiter, magdata, Weight, out muhat, out delta);
            for (int i = 0; i < M; i++)
                if (idx[i] == 0)
                    muhat[i] = (magdatatmp[i]) - (2 / magdatatmp[i]);
            for (int i = 0; i < M; i++)
                if (muhat[i] < 1e-7)
                    muhat[i] = 0;
            muhat = ILMath.multiplyElem(ILMath.sign((ILArray<double>)data), muhat).ToArray();
            double[] hugeMuInds = new double[M];
            for (int i = 0; i < M; i++)
                if (Math.Abs(muhat[i]) > Math.Abs(data[i]))
                    hugeMuInds[i] = 1;
                else
                    hugeMuInds[i] = 0;
            for (int i = 0; i < M; i++)
                if (hugeMuInds[i] == 1)
                    muhat[i] = data[i];


            return muhat;
        }

        private void btnInitialize_Click(object sender, RoutedEventArgs e)
        {
            lstInitialize.Items.Clear();
            BackgroundWorker worker = new BackgroundWorker();
            worker.WorkerReportsProgress = true;
            worker.DoWork += worker_InitializeDoWork;
            //worker.ProgressChanged += worker_ProgressChanged;

            worker.RunWorkerAsync();


        }
        void worker_InitializeDoWork(object sender, DoWorkEventArgs e)
        {
            //lstInitialize.Items.Clear();
            int progressPercentage = 0;
            R_sig = new double[] {12332276661.3281,
9947822939.95849,
7563369218.58887,
5178915497.21925,
30621245506.9055,
3120621622.01768,
25915965177.9550,
28619015346.9366,
110356303.151961,
22541284476.0305,
9575308198.13062,
8160683063.72067,
15967090494.4184,
2064971447.97539,
7848496134.88200,
5414271194.58297,
13449637961.5814,
12728850175.9023,
7645244387.08958,
19467429194.7894,
4296542716.76360,
12331704793.4057,
28350130551.1518,
1092207574.26882,
20815802108.8900,
19384069832.2560,
8904762656.46930,
13206844055.6089,
8574746213.94400,
7072983271.29331,
14531764554.1883,
4838479095.31312,
5087427804.14591,
1811709127.92492,
7307879492.93575,
12138006122.5229,
805656684.280176,
2146088831.76499,
4894287229.97550,
4910922111.41631,
403877998.166368,
364076746.813879,
414984886.549544,
950896493.292471,
566497536.201987,
1963938245.44071,
749507746.697584,
2116218177.08245,
1675908857.59278,
935305754.518634,
273474576.940248,
330377735.590040,
601781440.387666,
378260826.225747,
2616586210.41058,
181476148.251766,
629961724.333347,
1802887740.40988,
3517966021.12537,
654629016.878395,
3837051003.35208,
4560047806.76529,
1235899644.96884,
6806371897.71312,
8458848703.07143,
2397793193.63648,
3751218572.06078,
96038548.6139259,
7200786986.50766,
11823013550.0401,
2219603442.70747,
4190198217.49479,
8908372558.71633,
3990467371.85232,
2945375399.30515,
3043593778.87402,
5594419527.10096,
2446133465.47028,
2319166606.05123,
5647845622.46787,
2833477524.94254,
877102914.777345,
1954425017.69592,
900935971.540272,
483020111.044956,
1882913018.11158,
345632492.842714,
198916622.236928,
1904087784.99941,
1241021500.94707,
42870128.0105238,
1372567802.11437,
588038755.496795,
3698074876.11944,
583136031.899085,
1542170739.67918,
1625165035.50514,
8349055440.95197,
531957580.142910,
872197649.332370,
912383505.913620,
2019317361.26335,
500742667.044799,
1272535024.17488,
2278764259.80488,
2037026448.78911,
2233063227.60755,
152387796.586678,
2172753915.19694,
4544015644.52167,
1604933795.31650,
1765213744.03869,
3619301896.81528,
974330797.458942,
293786192.244344,
826776355.826075,
1596761522.91419,
1998337069.08333,
643439237.895726,
750623510.866383,
1017953500.57968,
477250804.725105,
1126776950.41138,
747120183.198999,
179593474.786727,
338189830.879543,
386049773.151436,
619280376.573806,
2588036234.84450,
5872823255.10722,
4961857655.29169,
106594008.009253,
934144123.677331,
3374921800.11116,
1267252058.29333,
406807795.754344,
369225273.401397,
1180042749.59626,
3643146880.52589,
2384883207.12815,
5003339283.36509,
1999348945.31098,
166083679.426466,
1477739310.93068,
1522372549.90211,
1802375216.24546,
1183376672.29359,
1546690576.91711,
74991760.5966050,
2040431390.58927,
2623557157.62268,
1815462631.81004,
3128180545.30722,
2175510731.29509,
1720833359.64052,
3269715300.90010,
2606325736.60032,
550572163.478881,
2264794728.49204,
1614966387.14112,
1891761993.32461,
2605104882.70652,
168970427.768299,
469307369.147672,
2508182882.92246,
828545630.084678,
1872172096.78369,
1853986893.78870,
1662936390.76838,
3675555064.06574,
4121954595.16545,
1058726053.53832,
2368190712.01297,
3775996587.28453,
2280732356.14948,
2270246971.92113,
578636714.709669,
4078021955.70376,
6375593831.57741,
739421695.621359,
5539245921.16653,
2161236346.32221,
10901237881.5075,
16750345873.2885,
3143667728.73852,
8213709982.68855,
9338761049.71514,
5331303829.98344,
9666296335.96934,
4855959443.44521,
5955990440.85444,
9709289288.79463,
3503845077.14861,
8431348634.82740,
6007091912.46386,
4570394882.36061,
9788211746.46656,
3941241174.24686,
2365274263.89020,
557079182.023114,
1075079957.71237,
3688126524.07716,
1519822928.35691,
196297701.041429,
4145784030.58913,
3866931037.29152,
347026150.166799,
2720024684.60798,
968194744.492103,
2649142081.91298,
7777708176.89416,
540641600.518077,
4360089943.87612,
4244116313.52459,
7913189598.41322,
7067455359.24602,
2187455245.05149,
16701000314.3405,
18533788212.3797,
20279088103.9418,
25548758888.3928,
13236852661.3248,
8505255038.09202,
13584513205.3807,
8654438076.15592,
12597556194.3331,
22612375584.7234,
30778392521.1511,
22302553148.7575,
24898287337.8017,
25483739699.6527,
5681876843.28218,
18969580042.5186,
23025510301.3854,
21965285520.5315,
35952374830.7482,
16354473780.6436,
12731052665.5670,
15684711704.0794,
12565916438.6092,
16244174384.0787,
7570466186.31500,
32876879062.1670,
16009227195.7853,
20419912862.4746,
53354283343.8735,
12969294591.0909,
21592205803.2188,
42723523753.7958,
2405429897.34576,
23449255461.0600,
16363942896.2186,
4950818040.12506,
19485637825.4370,
17101184104.0674,
14716730382.6977,
 };
            
            int len = R_sig.Length;
            for (int i = 0; i < len; i++)
            {
                progressPercentage = Convert.ToInt32(((double)(i + 1) / len) * 100);
                Dispatcher.Invoke(new Action(delegate ()
                {
                    lstInitialize.Items.Add(R_sig[i].ToString());
                }), System.Windows.Threading.DispatcherPriority.Normal);
                (sender as BackgroundWorker).ReportProgress(i);
                System.Threading.Thread.Sleep(10);
            }
        }

        void worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            //progressBar.Value = e.ProgressPercentage;
        }
        private void btnExecute_Click(object sender, RoutedEventArgs e)
        {
            lstListofModules.Items.Add("Algorithm Step 8. Optionally Remove 3/5 points spectral DC → R_sig Completed");
            //Algorithm Step 9. 	Initialize coefficients with modified 〖Sym〗_8 coefficients for a 5  level filtering → 〖nlvl〗_5
            String WavelettName = "Modified Sym 8"; int nlvl = 5;
            lstListofModules.Items.Add("Algorithm Step 9. Initialize coefficients with modified 〖Sym〗_8 coefficients for a 5  level filtering → 〖nlvl〗_5) Completed");

            //Algorithm Step 10. Calculate the Orthogonal wavelet filter set.
            double[] OrthloFit = null;
            if (WavelettName == "Modified Sym 8")
            {
                OrthloFit = new double[16] { 0.00133639669640, -0.00021419715012, -0.01057284326418, 0.00269319437688, 0.03474523295559, -0.01924676063167, -0.03673125438038, 0.25769933518654, 0.54955331526901, 0.34037267359439, -0.04332680770282, -0.10132432764282, 0.00537930587524, 0.02241181152181, -0.00038334544811, -0.00239172925575 };
                double sumorthloFit = OrthloFit.Sum();
                for (int i = 0; i < OrthloFit.Length; i++)
                    OrthloFit[i] = OrthloFit[i] / sumorthloFit;
            }
            lstListofModules.Items.Add("Algorithm Step 10. Calculate the Orthogonal wavelet filter set. completed");

            //Algorithm Step 10.1. Calculate Decomposition Low - pass, high - pass filter values → LoD,HiD
            //Algorithm Step 10.2. Calculate the Reconstruction Low-pass, high - pass filter values → LoR,HiR
            double[] LoD = new double[OrthloFit.Length];
            double[] HiD = new double[OrthloFit.Length];
            double[] LoR = new double[OrthloFit.Length];
            double[] HiR = new double[OrthloFit.Length];
            for (int i = 0; i < OrthloFit.Length; i++)
                LoR[i] = Math.Sqrt(2) * OrthloFit[i];

            for (int i = 0; i < LoR.Length; i++)
            {
                Double y = Convert.ToDouble(LoR[(LoR.Length - 1) - i]);
                if (i % 2 == 0) HiR[i] = y;
                else HiR[i] = -y;
            }

            for (int i = 0; i < HiR.Length; i++)
                HiD[i] = Convert.ToDouble(HiR[(HiR.Length - 1) - i]);
            for (int i = 0; i < LoR.Length; i++)
                LoD[i] = Convert.ToDouble(LoR[(LoR.Length - 1) - i]);
            lstListofModules.Items.Add("Algorithm Step 10.1. Calculate Decomposition Low - pass, high - pass filter values → LoD,HiD completed");
            lstListofModules.Items.Add("Algorithm Step 10.2. Calculate the Reconstruction Low-pass, high - pass filter values → LoR,HiR completed");

            //Algorithm Step 11.Calculate using eBayes taking Segment input:total Process
            lstListofModules.Items.Add("Algorithm Step 11. Calculate using eBayes taking Segment input:total Process completed");

            Double[] SegmentOutPut = MultiSignalWaveletDecomposition_eBayes(R_sig, nlvl, LoD, HiD, LoR, HiR);

            int len = SegmentOutPut.Length;
            for (int i = 0; i < len; i++)
            {
                lstExecute.Items.Add(SegmentOutPut[i].ToString());
            }
        }
        
    }
}
