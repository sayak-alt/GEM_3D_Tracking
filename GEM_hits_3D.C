#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TAxis.h"

void GEM_hits_3D() {
    const int nPoints = 1000000; // Number of points in the event

    // Create TGraph2D object
    TGraph2D *g3 = new TGraph2D();
    vector<double> x1, x2, x3;
    vector<double> y1, y2, y3;
    vector<double> z1, z2, z3;
    double a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6,b7;
    ifstream infile, infile1, infile2;
    //infile.open("corrected_x1x2x3.txt");
    infile.open("/Users/sayakchatterjee/Desktop/Sam_ug_2025/scripts/sorted_data_file_run_number_x1_x2_x3.txt");
    while(infile)
    {
        infile>>a1>>a2>>a3>>a4>>a5>>a6>>a7;
        if (a2>0 && a4>0 && a6>0 && a2<100 && a4<100 && a6<100)
        {
            x3.push_back(a2);
            z3.push_back(430);
            x2.push_back(a4);
            z2.push_back(230);
            x1.push_back(a6);
            z1.push_back(0);
        }
    }
    infile1.open("/Users/sayakchatterjee/Desktop/Sam_ug_2025/scripts/sorted_data_file_run_number_y1_y2_y3.txt");
    while(infile1)
    {
        infile1>>b1>>b2>>b3>>b4>>b5>>b6>>b7;
        if (b2>0 && b4>0 && b6>0 && b2<200 && b4<200 && b6<200)
        {
            y3.push_back(b2);
            y2.push_back(b4);
            y1.push_back(b6);
        }

    }


    // Example: fill with random (x, y, z) data representing hits
    for (int i = 0; i < x1.size(); ++i) {
        //double x = gRandom->Uniform(-10, 10);
       // double y = gRandom->Uniform(-10, 10);
        //double z = gRandom->Uniform(-10, 10);
        g3->SetPoint(i, x1[i], y1[i], z1[i]);
    }
    for (int i = 0; i < x2.size(); ++i) {
        //double x = gRandom->Uniform(-10, 10);
       // double y = gRandom->Uniform(-10, 10);
        //double z = gRandom->Uniform(-10, 10);
        g3->SetPoint(i+x1.size() , x2[i], y2[i], z2[i]);
    }
    for (int i = 0; i < x3.size(); ++i) {
        //double x = gRandom->Uniform(-10, 10);
       // double y = gRandom->Uniform(-10, 10);
        //double z = gRandom->Uniform(-10, 10);
        g3->SetPoint(i+x1.size()+x2.size(), x3[i], y3[i], z3[i]);
    }

    // Create a canvas to draw
    TCanvas *c1 = new TCanvas();
    g3->SetTitle("3D Event Hit Distribution;X [mm];Y [mm];Z [mm]");
    g3->GetXaxis()->SetRangeUser(0,100);
    g3->GetYaxis()->SetRangeUser(0,200);
    g3->GetZaxis()->SetRangeUser(0,700);

    g3->SetMarkerStyle(8);
    g3->SetMarkerColor(kBlue);
    g3->SetMarkerSize(0.2);
    g3->GetZaxis()->SetLimits(0,600);
    g3->Draw("P"); // P0 = points only, no interpolation

    // Keep canvas open
    c1->Update();
}
