#include <iostream>
#include <vector>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

void GEMs_residuals_x() {
    std::vector<float> z1, z2, z3;
    std::vector<float> x1, x2, x3;

    // Histograms
    TH1D *x1_uncorrected_angle_diff = new TH1D("x1_uncorrected_angle_diff", "X_{1} Difference", 200, -20, 20);
    TH1D *x2_uncorrected_angle_diff = new TH1D("x2_uncorrected_angle_diff", "X_{2} Difference", 200, -20, 20);
    TH1D *x3_uncorrected_angle_diff = new TH1D("x3_uncorrected_angle_diff", "X_{3} Difference", 200, -20, 20);
    
    TH1D *x1_corrected_angle_diff = new TH1D("x1_corrected_angle_diff", "X_{1} Difference", 200, -5, 5);
    TH1D *x2_corrected_angle_diff = new TH1D("x2_corrected_angle_diff", "X_{2} Difference", 200, -5, 5);
    TH1D *x3_corrected_angle_diff = new TH1D("x3_corrected_angle_diff", "X_{3} Difference", 200, -5, 5);

    // Reusable objects
    TH2D *histo_x1x2z1z2 = new TH2D("histo_x1x2z1z2", "X1/X2 vs Z1/Z2", 512, -0.5, 511.5, 1100, -0.5, 1099.5);
    TH2D *histo_x1x3z1z3 = new TH2D("histo_x1x3z1z3", "X1/X3 vs Z1/Z3", 512, -0.5, 511.5, 1100, -0.5, 1099.5);
    TH2D *histo_x2x3z2z3 = new TH2D("histo_x2x3z2z3", "X2/X3 vs Z2/Z3", 512, -0.5, 511.5, 1100, -0.5, 1099.5);

    TF1 *f1 = new TF1("f1", "pol1", 0, 512);
    f1->SetParameters(1, 1);
    TF1 *f2 = new TF1("f2", "pol1", 0, 512);
    f2->SetParameters(1, 1);
    TF1 *f3 = new TF1("f3", "pol1", 0, 512);
    f3->SetParameters(1, 1);
    
    TF1 *f11 = new TF1("f1", "gaus", -15, 15);
    f1->SetParameters(1, 1);
    TF1 *f21 = new TF1("f2", "gaus", -15, 15);
    f2->SetParameters(1, 1);
    TF1 *f31 = new TF1("f3", "gaus", -15, 15);
    f3->SetParameters(1, 1);

    // Read data from file
    std::ifstream infile("sorted_data_file_run_number_x1_x2_x3.txt");
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
        return;
    }

    float a1, a2, a3, a4, a5, a6, a7;
    while (infile >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7) {
        x1.push_back(a2);
        //z1.push_back(a3);
        z1.push_back(430);
        x2.push_back(a4);
        //z2.push_back(a5);
        z2.push_back(230);
        x3.push_back(a6+12.14);
        //z3.push_back(a7);
        z3.push_back(0);
    }
    infile.close();

    // Processing loop
    const size_t N = x1.size();
    for (size_t k = 0; k < N; ++k) {
    //  for (size_t k = 0; k < 1000; ++k) {
        // First fit with x1, x2
        histo_x1x3z1z3->Reset();
        histo_x1x3z1z3->Fill(x1[k] , z1[k]);
        histo_x1x3z1z3->Fill(x3[k] , z3[k]);
        histo_x1x3z1z3->Fit(f2, "QRN");
        
        histo_x1x2z1z2->Reset();
        histo_x1x2z1z2->Fill(x1[k] , z1[k]);
        histo_x1x2z1z2->Fill(x2[k] , z2[k]);
        histo_x1x2z1z2->Fit(f3, "QRN");
        
        histo_x2x3z2z3->Reset();
        histo_x2x3z2z3->Fill(x2[k] , z2[k]);
        histo_x2x3z2z3->Fill(x3[k] , z3[k]);
        histo_x2x3z2z3->Fit(f1, "QRN");
        cout << " event no " << k << " total no of events " << N << endl;
        
        double x2_cal = (230.0 - f2->GetParameter(0)) / f2->GetParameter(1);
        double x1_cal = (430.0 - f1->GetParameter(0)) / f1->GetParameter(1);
        double x3_cal = (0.0 - f3->GetParameter(0)) / f3->GetParameter(1);
        
        //....X....///
        // double x3_diff = x2[k] - x3_cal -4.5;
        // double x3_diff = x3[k] - x3_cal ;
        // double x3_diff = x1[k] - x3_cal ;
        //...Y.....//
        //   double x3_diff = x1[k] - x3_cal + 0.23 ;
        double x1_diff = x1[k] - x1_cal;
        double x2_diff = x2[k] - x2_cal;
        double x3_diff = x3[k] - x3_cal;
        
        x1_uncorrected_angle_diff->Fill(x1_diff);
        x1_uncorrected_angle_diff->Fit(f11,"QRN");
        double x1_offset = f11->GetParameter(1);
        double x1_diff_cor = x1_diff - x1_offset;
        x1_corrected_angle_diff->Fill(x1_diff_cor);
        
        
        x2_uncorrected_angle_diff->Fill(x2_diff);
        x2_uncorrected_angle_diff->Fit(f21,"QRN");
        double x2_offset = f21->GetParameter(1);
        double x2_diff_cor = x2_diff - x2_offset;
        x2_corrected_angle_diff->Fill(x2_diff_cor);
        
        
        
        
        
        
        x3_uncorrected_angle_diff->Fill(x3_diff);
        x3_uncorrected_angle_diff->Fit(f31,"QRN");
        double x3_offset = f31->GetParameter(1);
        double x3_diff_cor = x3_diff - x3_offset;
        x3_corrected_angle_diff->Fill(x3_diff_cor);
        
    }
    
    x1_corrected_angle_diff->Fit(f11,"R");
    cout << " Mean for corrected X1 " << f11->GetParameter(0) << " Sigma " << f11->GetParameter(0) << endl;
    x2_corrected_angle_diff->Fit(f21,"R");
    cout << " Mean for corrected X2 " << f21->GetParameter(0) << " Sigma " << f21->GetParameter(0) << endl;
    x3_corrected_angle_diff->Fit(f31,"R");
    cout << " Mean for corrected X3 " << f31->GetParameter(0) << " Sigma " << f31->GetParameter(0) << endl;

    // Plot results
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->Divide(3,1);
    x3_uncorrected_angle_diff->GetXaxis()->SetTitle("X_{3}^{measured} - X_{3}^{calculated} (mm)");
    x3_uncorrected_angle_diff->GetYaxis()->SetTitle("Counts");
    x1_uncorrected_angle_diff->GetXaxis()->SetTitle("X_{1}^{measured} - X_{1}^{calculated} (mm)");
    x1_uncorrected_angle_diff->GetYaxis()->SetTitle("Counts");
    x2_uncorrected_angle_diff->GetXaxis()->SetTitle("X_{2}^{measured} - X_{2}^{calculated} (mm)");
    x2_uncorrected_angle_diff->GetYaxis()->SetTitle("Counts");
    c1->cd(1);
    x1_uncorrected_angle_diff->Draw();
    c1->cd(2);
    x2_uncorrected_angle_diff->Draw();
    c1->cd(3);
    x3_uncorrected_angle_diff->Draw();
    
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->Divide(3,1);
    x3_corrected_angle_diff->GetXaxis()->SetTitle("X_{3}^{measured} - X_{3}^{calculated} (mm)");
    x3_corrected_angle_diff->GetYaxis()->SetTitle("Counts");
    x1_corrected_angle_diff->GetXaxis()->SetTitle("X_{1}^{measured} - X_{1}^{calculated} (mm)");
    x1_corrected_angle_diff->GetYaxis()->SetTitle("Counts");
    x2_corrected_angle_diff->GetXaxis()->SetTitle("X_{2}^{measured} - X_{2}^{calculated} (mm)");
    x2_corrected_angle_diff->GetYaxis()->SetTitle("Counts");
    c2->cd(1);
    x1_corrected_angle_diff->Draw();
    c2->cd(2);
    x2_corrected_angle_diff->Draw();
    c2->cd(3);
    x3_corrected_angle_diff->Draw();


}

