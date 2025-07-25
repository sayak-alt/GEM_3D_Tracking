// This is the part of the code which does the independent 2D fits into a 3D fit
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <ctime>

#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TLegend.h"
#include "TFile.h"
#include "TView.h"
#include "TText.h"

struct TrackEvent {
    int run;
    double x[3], y[3], z[3];
    double angleDeg;
    double chi2_ndf;
};

void fit_and_angle_from_points(double* x, double* y, double* z, double& dx_dz, double& dy_dz, double& chi2, double& signed_angle_deg) {
    TGraphErrors gx(3, z, x, 0, 0);
    TGraphErrors gy(3, z, y, 0, 0);

    TF1 fx("fx", "pol1");
    TF1 fy("fy", "pol1");
    gx.Fit(&fx, "Q0");
    gy.Fit(&fy, "Q0");

    dx_dz = fx.GetParameter(1);
    dy_dz = fy.GetParameter(1);

    double chi2x = fx.GetChisquare();
    double chi2y = fy.GetChisquare();
    chi2 = (chi2x + chi2y) / 2.0;

    // Signed polar angle in x-z plane (using slope dx/dz)
    signed_angle_deg = std::atan(dx_dz) * 180.0 / TMath::Pi();
}

void muon_tracker() {
    std::ifstream fx("sorted_data_file_run_number_x1_x2_x3.txt");
    std::ifstream fy("sorted_data_file_run_number_y1_y2_y3.txt");

    if (!fx.is_open() || !fy.is_open()) {
        std::cerr << "Error opening input files!" << std::endl;
        return;
    }

    std::map<int, TrackEvent> events;
    const double zGEM[3] = {430.0, 230.0, 0.0};

    int run;
    double x1, z1x, x2, z2x, x3, z3x;
    while (fx >> run >> x1 >> z1x >> x2 >> z2x >> x3 >> z3x) {
        TrackEvent evt;
        evt.run = run;
        evt.x[0] = x1;
        evt.x[1] = x2;
        evt.x[2] = x3 ;
        evt.z[0] = zGEM[0];
        evt.z[1] = zGEM[1];
        evt.z[2] = zGEM[2];
        events[run] = evt;
    }

    double y1, z1y, y2, z2y, y3, z3y;
    while (fy >> run >> y1 >> z1y >> y2 >> z2y >> y3 >> z3y) {
        if (events.find(run) != events.end()) {
            events[run].y[0] = y1;
            events[run].y[1] = y2;
            events[run].y[2] = y3 ;
        }
    }

    // Histogram for signed angle in x-z plane
    TH1D* hSignedAngleXZ = new TH1D("hSignedAngleXZ", "Signed Angle in X-Z Plane;Angle [deg];Counts", 360, -90, 90);

    // Fit and calculate signed angle for all events
    std::vector<int> keys;
    for (auto& [run, evt] : events) {
        double dx_dz, dy_dz, chi2, signed_angle;
        fit_and_angle_from_points(evt.x, evt.y, evt.z, dx_dz, dy_dz, chi2, signed_angle);
        evt.angleDeg = signed_angle;
        evt.chi2_ndf = chi2;
        hSignedAngleXZ->Fill(signed_angle);
        keys.push_back(run);
    }

    // Randomly sample 10 events for 3D plotting
    std::shuffle(keys.begin(), keys.end(), std::default_random_engine(static_cast<unsigned>(time(nullptr))));
    TCanvas* c3D = new TCanvas("c3D", "Random 10 Tracks in 3D", 1000, 800);
    TView* view = TView::CreateView(1);
    view->SetRange(0, -300, -300, 450, 300, 300);
    view->SetPerspective();
    view->ShowAxis();

    TLegend* legend = new TLegend(0.75, 0.80, 0.95, 0.92);
    legend->AddEntry((TObject*)0, "Red: Hits", "");
    legend->AddEntry((TObject*)0, "Blue: Fit", "");

    for (int i = 0; i < 10 && i < (int)keys.size(); ++i) {
        auto& evt = events[keys[i]];
        double dx_dz, dy_dz, chi2, signed_angle;
        fit_and_angle_from_points(evt.x, evt.y, evt.z, dx_dz, dy_dz, chi2, signed_angle);

        // Points
        TPolyMarker3D* points = new TPolyMarker3D(3);
        for (int j = 0; j < 3; ++j)
            points->SetPoint(j, evt.x[j], evt.y[j], evt.z[j]);
        points->SetMarkerStyle(20);
        points->SetMarkerColor(kRed);
        points->SetMarkerSize(1.0);
        points->Draw();

        // Fit line
        TPolyLine3D* line = new TPolyLine3D(2);
        for (int j = 0; j < 2; ++j) {
            double z = (j == 0) ? 0.0 : 450.0;
            double x = dx_dz * z + (evt.x[2] - dx_dz * evt.z[2]);
            double y = dy_dz * z + (evt.y[2] - dy_dz * evt.z[2]);
            line->SetPoint(j, x, y, z);
        }
        line->SetLineColor(kBlue);
        line->SetLineWidth(1);
        line->Draw("same");
    }
    legend->Draw();

    // Draw signed angle histogram
    TCanvas* cAngle = new TCanvas("cAngle", "Signed Angle X-Z Histogram", 800, 600);
    hSignedAngleXZ->SetLineColor(kBlue + 2);
    hSignedAngleXZ->SetLineWidth(2);
    hSignedAngleXZ->Draw();
}


