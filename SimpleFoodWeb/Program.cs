//
// Simple food chain/food web model
//
// Original code by Mick Follows (14/12/15). Described in Knoll & Follows, Royal Proc B (2016).
// Simpler version of model described in Ward et al., L&O (2013), Ward et al., JPR (2014) and Ward & Follows (2016).
// All after Armstrong (1994), see Ward & Follows, JPR (2014).
//
// Ported from  python 2.7 to C# by Erik A. Mousing (13/11/16)
//


using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace SimpleFoodWeb
{
    class Program
    {
        static void Main(string[] args)
        {
            // set number of trophic levels
            var nLevels = 8;

            // Flag for Holling type I or Holling type II grazing 
            // 1 for Holling I, 2 for Holling II
            // ONLY HOLLING TYPE I GRAZING HAS BEEN IMPLEMENTED
            var holling = 1; // Do not change this parameter

            // set time step parameters
            var dt = 0.05;
            var maxTime = 3000.0;
            var nStepMax = Convert.ToInt32(Math.Round(maxTime / dt));
            var timeOut = 1000.0;
            var nStepOut = nStepMax * timeOut / maxTime;
            var nStepOutMax = Convert.ToInt32(nStepMax / nStepOut);

            // Initialize working arrays
            var bmassN = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            var dBmassNDt = new double[nLevels];

            // Initialize output arrays
            var bmassN_Out = new double[nLevels, nStepOutMax];
            var autotrophy_Out = new double[nLevels, nStepOutMax];
            var heterotrophy_Out = new double[nLevels, nStepOutMax];
            var predation_Out = new double[nLevels, nStepOutMax];
            var respiration_Out = new double[nLevels, nStepOutMax];
            var time_Out = new double[nStepOutMax];
            var nitrate_Out = new double[nStepOutMax];

            // Define cell volumes
            var cellVol = new double[nLevels];
            cellVol[0] = Math.Pow(10.0, 1.0);
            cellVol[1] = Math.Pow(10.0, 1.0);
            cellVol[2] = Math.Pow(10.0, 1.5);
            cellVol[3] = Math.Pow(10.0, 1.5);
            cellVol[4] = Math.Pow(10.0, 2.0);
            cellVol[5] = Math.Pow(10.0, 2.0);
            cellVol[6] = Math.Pow(10.0, 2.5);
            cellVol[7] = Math.Pow(10.0, 2.5);

            // Print cell volumes to screen
            Console.WriteLine("[{0}]", string.Join(", ", cellVol));

            // Insert code that prints cell volumes to a file

            // Initialize arrays for cell quotas
            var quotaCarbon = new double[nLevels];
            var quotaNitrate = new double[nLevels];

            var aCell = 18.7;
            var bCell = 0.89;

            // Calculate carbon cell quota for each size group
            for (var i = 0; i < quotaCarbon.Length; i++)
            {
                quotaCarbon[i] = aCell * Math.Pow(cellVol[i], bCell);
            }

            // Calculate nitrate cell quota for each size group based on Redfield proportions
            for (var i = 0; i < quotaNitrate.Length; i++)
            {
                quotaNitrate[i] = quotaCarbon[i] * (16.0 / 106.0);
            }

            // Convert nitrate quota to micromoles N cell-1
            for (var i = 0; i < quotaNitrate.Length; i++)
            {
                quotaNitrate[i] = quotaNitrate[i] * 1.0e6 / 1.0e15;
            }

            // Initiate arrays for autotrophic traits
            var vmaxN = new double[nLevels];
            var kn = new double[nLevels];
            var kResp = new double[nLevels];
            var specVmaxN = new double[nLevels];

            // Set even numbered types to be obligate autotrophs
            vmaxN[0] = 9.1e-9 * Math.Pow(cellVol[0], 0.67);
            vmaxN[2] = 9.1e-9 * Math.Pow(cellVol[2], 0.67);
            vmaxN[4] = 9.1e-9 * Math.Pow(cellVol[4], 0.67);
            vmaxN[6] = 9.1e-9 * Math.Pow(cellVol[6], 0.67);

            // Set odd numbered types to be obligate heterotrophs
            vmaxN[1] = 0.0;
            vmaxN[3] = 0.0;
            vmaxN[5] = 0.0;
            vmaxN[7] = 0.0;

            // Set half-saturation constant for all types
            for (var i = 0; i < kn.Length; i++)
            {
                kn[i] = 0.17 * Math.Pow(cellVol[i], 0.27);
            }

            // Set mortality term to be the same for all types
            for (var i = 0; i < kResp.Length; i++)
            {
                kResp[i] = 0.03;
            }

            //
            for (var i = 0; i < specVmaxN.Length; i++)
            {
                specVmaxN[i] = vmaxN[i] / quotaNitrate[i];
            }

            // Print autotrophic trait arrays to console
            Console.WriteLine("vmaxN:");
            Console.WriteLine("[{0}]", string.Join(", ", vmaxN));
            Console.WriteLine("specVmaxN:");
            Console.WriteLine("[{0}]", string.Join(", ", specVmaxN));
            Console.WriteLine("kn:");
            Console.WriteLine("[{0}]", string.Join(", ", kn));

            // Initiate arrays for heterotrophic traits
            var gmax = new double[nLevels, nLevels];
            var gmax1 = new double[nLevels, nLevels];
            var kbN = new double[nLevels, nLevels];
            var gamma = new double[nLevels, nLevels];
            var eye = new double[nLevels, nLevels];
            var grazInteract = new double[nLevels, nLevels];

            // Populate grazing interaction matrix

            grazInteract[0, 1] = 1.0;
            grazInteract[2, 3] = 1.0;
            grazInteract[4, 5] = 1.0;
            grazInteract[6, 7] = 1.0;

            // Figure out how to print a matrix to the console

            var aGraz = 0.5;
            var bGraz = -0.16;

            // Calculate maximum grazing rate for all types
            for (var i = 0; i < gmax.GetLength(0); i++)
            {
                for (var j = 0; j < gmax.GetLength(1); j++)
                {
                    gmax[i, j] = aGraz * Math.Pow(cellVol[j], bGraz) * grazInteract[i, j];
                    eye[i, j] = i;
                }
            }

            // Set gmax1 = grazing eate for Holling I ((mol N l-1 day)-1)
            for (var i = 0; i < gmax1.GetLength(0); i++)
            {
                for (var j = 0; j < gmax1.GetLength(1); j++)
                {
                    // Ask Mick about this.. gmax and gmax1 are identical(?)
                    gmax1[i, j] = gmax[i, j] / 1.0;
                }
            }

            // Print matrix to console
            //
            //var rowCount = gmax.GetLength(0);
            //var colCount = gmax.GetLength(1);
            //for (int row = 0; row < rowCount; row++)
            //{
            //    for (int col = 0; col < colCount; col++)
            //        Console.Write(String.Format("{0}\t", gmax[row, col]));
            //    Console.WriteLine();
            //}

            // Calculate half-saturation constant for grazing (from Ward et al. (2013))
            for (var i = 0; i < kbN.GetLength(0); i++)
            {
                for (var j = 0; j < kbN.GetLength(1); j++)
                {
                    kbN[i, j] = 0.5 * (16.0 / 106.0) * grazInteract[i, j]; // micromol N biomass L-1 
                }
            }

            // Set trophic transfer efficiency
            for (var i = 0; i < gamma.GetLength(0); i++)
            {
                for (var j = 0; j < gamma.GetLength(1); j++)
                {
                    gamma[i, j] = 0.1 * grazInteract[i, j];
                }
            }

            // Working arrays and variables
            var heterotrophy = new double[nLevels];
            var autotrophy = new double[nLevels];
            var predation = new double[nLevels];
            var respiration = new double[nLevels];
            var imax = nLevels;
            var jmax = nLevels;
            var kmax = nLevels;

            // Start loop over snitrate (nitrate supply rate)

            for (var sStep = 0; sStep < 10; sStep++)
            {
                var sNitrate = 0.02 + sStep * 0.06;
                Console.WriteLine("sNitrate = {0}", sNitrate); // Place holder is not working
                // Initiate time related variables within the loop
                var nCount = 0;
                var nOut = 0;
                var time = 0.0;
                var dBmassNDtOld = 0.0;
                var dNitrateDtOld = 0.0;
                var dNitrateDt = 0.0;
                var prod = new double[nLevels];

                // Re-initialize biomass (convert to micromoles N l-1??)
                for (var i = 0; i < bmassN.Length; i++)
                {
                    bmassN[i] = bmassN[i] * 1.0e-2;
                }

                // Set nitrate concentration
                var nitrate = 1.0;

                // Start main time loop
                for (var nStep = 0; nStep < nStepMax; nStep++)
                {
                    // Evaluate autotrophy and respiration
                    for (var i = 0; i < autotrophy.Length; i++)
                    {
                        autotrophy[i] = (vmaxN[i] / quotaNitrate[i]) * (nitrate / (nitrate + kn[i])) * bmassN[i];
                    }

                    // Maintain respiration/background loss
                    for (var i = 0; i < respiration.Length; i++)
                    {
                        respiration[i] = kResp[i] * bmassN[i];
                    }

                    // Evaluate heterotrophy
                    for (var i = 0; i < imax; i++)
                    {
                        if (holling == 1)
                        {
                            // 
                            var sum1 = 0.0;
                            for (var j = 0; j < jmax; j++)
                            {
                                sum1 = sum1 + (gamma[j, i] * gmax1[j, i] * bmassN[i] * bmassN[j]);
                                heterotrophy[i] = sum1;
                            }

                            // Predation by all others (ask Mick)
                            var sum2 = 0.0;
                            for (var k = 0; k < kmax; k++)
                            {
                                sum2 = sum2 + (gmax1[i, k] * bmassN[k] * bmassN[i]);
                                predation[i] = sum2;
                            }
                        }
                    }


                    // To avoid numerical problems, do not let biomass decline below a low threshold value
                    for (var i = 0; i < nLevels; i++)
                    {
                        if (bmassN[i] < 1.0e-25)
                        {
                            respiration[i] = 0.0;
                            predation[i] = 0.0;
                        }
                    }

                    // Evaluate rates of change
                    dNitrateDt = -autotrophy.Sum() + sNitrate;

                    for (var i = 0; i < dBmassNDt.Length; i++)
                    {
                        dBmassNDt[i] = autotrophy[i] + heterotrophy[i] - respiration[i] - predation[i];
                    }

                    // Store total growth
                    for (var i = 0; i < prod.Length; i++)
                    {
                        prod[i] = autotrophy[i] + heterotrophy[i];
                    }

                    // Euler forward step
                    var bmassNNew = new double[nLevels];
                    double nitrateNew;

                    for (var i = 0; i < bmassNNew.Length; i++)
                    {
                        bmassNNew[i] = bmassN[i] + dBmassNDt[i] * dt;
                        bmassN[i] = bmassNNew[i];
                    }

                    nitrateNew = nitrate + dNitrateDt * dt;
                    nitrate = nitrateNew;

                    // Increment time
                    time = time + dt;

                    // Print to screen and put diagnostics into an array
                    nCount = nCount + 1;
                    if (nCount == nStepOut)
                    {
                        // Print values to screen
                        string screenOut1 = "";
                        screenOut1 += "Time: " + Math.Round(time, 0)
                            + "; Nitrate: " + nitrate.ToString("0.00E0")
                            + "; Biomass grp0:" + bmassN[0].ToString("0.00E0")
                            + "; Biomass grp1:" + bmassN[1].ToString("0.00E0")
                            + "; Biomass grp2:" + bmassN[2].ToString("0.00E0")
                            + "; Biomass grp3:" + bmassN[3].ToString("0.00E0");

                        Console.WriteLine(screenOut1);

                        // nCount is the counter to decide whether to write output - reset here
                        nCount = 0;

                        // nOut is the position in the output array
                        time_Out[nOut] = time;
                        nitrate_Out[nOut] = nitrate;
                        for (var i = 0; i < bmassN.Length; i++)
                        {
                            bmassN_Out[i, nOut] = bmassN[i];
                            autotrophy_Out[i, nOut] = autotrophy[i];
                            heterotrophy_Out[i, nOut] = heterotrophy[i];
                            predation_Out[i, nOut] = predation[i];
                            respiration_Out[i, nOut] = respiration[i];
                        }

                        // Increment the output counter
                        nOut = nOut + 1;
                    }
                }
                // End main time loop
                Console.WriteLine("end off loop");



                if (sStep == 0)
                {
                    Console.WriteLine("The first loop has now finished");

                }

                else
                {
                    Console.WriteLine("looping");
                }
            }
        }
    }
}
