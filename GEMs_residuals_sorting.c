#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

void GEMs_residuals_sorting() {
    int counter = 0;
    std::vector<float> z1_con, z2_con, z3_con;
    std::vector<float> x1_con, x2_con, x3_con;
    std::vector<float> run_number_sorted;
    
    double x1_corrected, x2_corrected, x3_corrected;
    int counter_gem1=0, counter_gem2=0, counter_gem3=0;
    
    std::vector<float> run_number_x1, run_number_x2, run_number_x3;
    std::vector<float> x1, x2, x3;
    
    std::ifstream infile1, infile2, infile3;
    
    // GEM-I
    infile1.open("output_file_run_3766_3768_x1y1_GEM-I.txt");
    float a1, a2, a3, a4, a5;
    while (infile1 >> a1 >> a2 >> a3 >> a4 >> a5)
    {
       // if (counter_gem1> 80000) break;
        run_number_x1.push_back(a1);
        x1.push_back(a2 * 0.39062500);
        counter_gem1++;

    }
    
    // GEM-II
    infile2.open("output_file_run_3766_3768_x2y2_GEM-II.txt");
    float b1, b2, b3, b4, b5;
    while (infile2 >> b1 >> b2 >> b3 >> b4 >> b5)
    {
        //if (counter_gem2 > 80000) break;
        run_number_x2.push_back(b1);
        x2.push_back(b2 * 0.39062500);
        counter_gem2++;

    }
    
    // GEM-III
    infile3.open("output_file_run_3766_3768_x3y3_GEM-III.txt");
    float r1, r2, r3, r4, r5;
    while (infile3 >> r1 >> r2 >> r3 >> r4 >> r5)
    {
       // if (counter_gem3 > 80000) break;
        run_number_x3.push_back(r1);
        x3.push_back(r2 * 0.39062500);
        counter_gem3++;
    }
    
    std::vector<int> run_number;
   // run_number.reserve(run_number_x1.size()); // Pre-allocate space for efficiency
   
    ofstream outfile;
    outfile.open("sorted_data_file_run_number_x1_x2_x3.txt");
    
    for (int j=0; j<run_number_x1.size(); j++)
    {
        for (int k=0; k<run_number_x2.size(); k++)
        {
            if (run_number_x2[k] == run_number_x1[j])
            {
                for (int l=0; l<run_number_x3.size(); l++)
                {
                    if (run_number_x3[l] == run_number_x1[j])
                    {
                        x1_con.push_back(x1[j]);
                        x2_con.push_back(x2[k]);
                        x3_con.push_back(x3[l]);
                        z3_con.push_back(0);          // Treating GEM-III z as the origin
                        z2_con.push_back(230);      // Distance between GEM-III & GEM-II
                        z1_con.push_back(430);      // Distance of GEM-I from GEM-III
                        std::cout << "run number " << run_number_x3[counter] << " x1 " << x1_con[counter] << " x2 " << x2_con[counter] << " x3 " << x3_con[counter] << std::endl;
                        outfile << run_number_x3[counter] << " \t " << x1_con[counter] <<  " \t " << z1_con[counter] << " \t " << x2_con[counter] << " \t " << z2_con[counter] << " \t "  << x3_con[counter] << " \t " << z3_con[counter] << endl;
                        counter++;
                    }
                }
            }
        }
    }
     
    outfile << " Run number " << " \t " << " x1 " << " \t " << " z1 " << " \t " << " x2 " << " \t " << " z2 " << " \t " << " x3 " << " \t " << " z3 " << endl;
    
    
}

