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

            // set time step parameters
            var dt = 0.02;
            var maxTime = 3000.0;
            var nStepMax = Convert.ToInt32(Math.Round(maxTime / dt));
            var timeOut = 1000.0;
            var nStepOut = nStepMax * timeOut / maxTime;
            var nStepOutMax = Convert.ToInt32(nStepMax / nStepOut);

            // Initialize working arrays
            var bmassN = new double[nLevels];
            var dBmassNdt = new double[nLevels];

            // Initialize output arrays
            var bmassNOut = new double[nLevels, nStepOutMax];
            var autotrophyOut = new double[nLevels, nStepOutMax];
            var heterotrophyOut = new double[nLevels, nStepOutMax];
            var predationOut = new double[nLevels, nStepOutMax];
            var respirationOut = new double[nLevels, nStepOutMax];
            var timOut = new double[nStepOutMax];
            var nitrateOut = new double[nStepOutMax];

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
            for(var i = 0; i < quotaNitrate.Length; i++)
            {
                quotaNitrate[i] = quotaCarbon[i] * (16.0 / 106.0);
            }

            // Convert nitrate quota to micromoles N cell-1
            for(var i = 0; i < quotaNitrate.Length; i++)
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
            for(var i = 0; i < kn.Length; i++)
            {
                kn[i] = 0.17 * Math.Pow(cellVol[i], 0.27);
            }

            // Set mortality term to be the same for all types
            for(var i = 0; i < kResp.Length; i++)
            {
                kResp[i] = 0.03;
            }

            //
            for(var i = 0; i < specVmaxN.Length; i++)
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
            for(var i = 0; i < gmax.GetLength(0); i++)
            {
                for(var j = 0; j < gmax.GetLength(1); j++)
                {
                    gmax[i, j] = aGraz * Math.Pow(cellVol[j], bGraz) * grazInteract[i, j];
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
            for(var i = 0; i < gamma.GetLength(0); i++)
            {
                for(var j = 0; j < gamma.GetLength(1); j++)
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
            
            for(var sStep = 0; sStep < 10; sStep++)
            {
                var sNitrate = 0.02 + sStep * 0.06;
                Console.WriteLine("sNitrate = {{1}}", sNitrate); // Place holder is not working
                // Initiate time related variables within the loop
                var nCount = 0;
                var nOut = 0;
                var time = 0.0;
                var dBmassNDtOld = 0.0;
                var dNitrateDtOld = 0.0;
                // Re-initialize biomass (convert to micromoles N l-1??)
                for(var i = 0; i < bmassN.Length; i++)
                {
                    bmassN[i] = bmassN[i] * 1.0e-2;
                }
                
                // Set nitrate concentration
                var nitrate = 1.0;

                // Start main time loop
                for(var nStep = 0; nStep < nStepMax; nStep++)
                {
                    // Evaluate autotrophy and respiration
                    for(var i = 0; i < autotrophy.Length; i++)
                    {
                        autotrophy[i] = (vmaxN[i] / quotaNitrate[i]) * (nitrate / (nitrate + kn[i])) * bmassN[i];
                    }

                    // Maintain respiration/background loss
                    for(var i = 0; i < respiration.Length; i++)
                    {
                        respiration[i] = kResp[i] * bmassN[i];
                    }

                    // Evaluate heterotrophy
                    // Only Holling I type has been implemented
                    //for(var i = 0; i < imax; i++)
                    //{
                    //    var sum1 = 0.0;
                    //    for(var j = 0; j < jmax; j++)
                    //    {
                            
                    //    }
                    //} 
                }
            }


        }
    }
}
