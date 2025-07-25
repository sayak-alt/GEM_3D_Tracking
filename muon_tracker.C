
// This is the part of the code which does the 3D tracking
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
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

struct TrackEvent {
    int run;
    double x[3], y[3], z[3];
    double angleDeg;
    double chi2_ndf;
};

void fit3DLine(const double* x, const double* y, const double* z, double centroid[3], double direction[3]) {
    const int N = 3;
    centroid[0] = centroid[1] = centroid[2] = 0.0;
    for (int i = 0; i < N; ++i) {
        centroid[0] += x[i];
        centroid[1] += y[i];
        centroid[2] += z[i];
    }
    for (int i = 0; i < 3; ++i)
        centroid[i] /= N;

    TMatrixD cov(3, 3);
    cov.Zero();
    for (int i = 0; i < N; ++i) {
        double dx = x[i] - centroid[0];
        double dy = y[i] - centroid[1];
        double dz = z[i] - centroid[2];
        cov(0,0) += dx*dx; cov(0,1) += dx*dy; cov(0,2) += dx*dz;
        cov(1,0) += dy*dx; cov(1,1) += dy*dy; cov(1,2) += dy*dz;
        cov(2,0) += dz*dx; cov(2,1) += dz*dy; cov(2,2) += dz*dz;
    }
    cov *= (1.0 / N);

    TMatrixDSym covSym(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            covSym(i, j) = cov(i, j);

    TMatrixDSymEigen eigenDecomp(covSym);
    TVectorD eigenValues = eigenDecomp.GetEigenValues();
    TMatrixD eigenVectors = eigenDecomp.GetEigenVectors();

    int maxIndex = 0;
    double maxValue = eigenValues[0];
    for (int i = 1; i < 3; ++i) {
        if (eigenValues[i] > maxValue) {
            maxValue = eigenValues[i];
            maxIndex = i;
        }
    }

    for (int i = 0; i < 3; ++i)
        direction[i] = eigenVectors(i, maxIndex);

    double mag = std::sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
    for (int i = 0; i < 3; ++i)
        direction[i] /= mag;
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
            events[run].y[2] = y3;
        }
    }

    TH1D* hPolarAngle = new TH1D("hPolarAngle", "Signed Polar Angle (3D Fit);Polar Angle [deg];Counts", 360, -90, 90);

    std::vector<int> keys;
    for (auto& [run, evt] : events) {
        double centroid[3], direction[3];
        fit3DLine(evt.x, evt.y, evt.z, centroid, direction);

        // Compute polar angle (signed angle in X-Z plane)
        double angle_rad = std::atan2(direction[0], direction[2]); // x/z
        double angle_deg = angle_rad * 180.0 / TMath::Pi();
        evt.angleDeg = angle_deg;
        hPolarAngle->Fill(angle_deg);
        keys.push_back(run);
    }

    // Randomly sample 10 events to draw in 3D
    std::shuffle(keys.begin(), keys.end(), std::default_random_engine(static_cast<unsigned>(time(nullptr))));
    TCanvas* c3D = new TCanvas("c3D", "Random 10 Tracks in 3D", 1000, 800);
    TView* view = TView::CreateView(1);
    view->SetRange(0, -300, -300, 450, 300, 300);
    view->SetPerspective();
    view->ShowAxis();

    TLegend* legend = new TLegend(0.75, 0.80, 0.95, 0.92);
    legend->AddEntry((TObject*)0, "Red: Hits", "");
    legend->AddEntry((TObject*)0, "Blue: 3D Fit", "");

    for (int i = 0; i < 10 && i < keys.size(); ++i) {
        auto& evt = events[keys[i]];
        double centroid[3], direction[3];
        fit3DLine(evt.x, evt.y, evt.z, centroid, direction);

        TPolyMarker3D* points = new TPolyMarker3D(3);
        for (int j = 0; j < 3; ++j)
            points->SetPoint(j, evt.x[j], evt.y[j], evt.z[j]);
        points->SetMarkerStyle(20);
        points->SetMarkerColor(kRed);
        points->SetMarkerSize(1.0);
        points->Draw();

        TPolyLine3D* line = new TPolyLine3D(2);
        for (int j = 0; j < 2; ++j) {
            double len = (j == 0) ? -250.0 : 250.0;
            double x = centroid[0] + len * direction[0];
            double y = centroid[1] + len * direction[1];
            double z = centroid[2] + len * direction[2];
            line->SetPoint(j, x, y, z);
        }
        line->SetLineColor(kBlue);
        line->SetLineWidth(1);
        line->Draw("same");
    }
    legend->Draw();

    // Plot polar angle histogram
    TCanvas* cAngle = new TCanvas("cAngle", "Polar Angle Histogram", 800, 600);
    hPolarAngle->SetLineColor(kBlue + 2);
    hPolarAngle->SetLineWidth(2);
    hPolarAngle->Draw();
}

