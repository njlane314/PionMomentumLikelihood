//root -l -b -q 'PlotValidation.c("../pimom_output.root")'
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

void PlotValidation(const char* filename) {
    TFile* file = new TFile(filename);
    TTree* tree = (TTree*)file->Get("ana/Tree");

    double sim_px, sim_py, sim_pz, sim_purity, sim_completeness;
    std::vector<double> *sp_x = nullptr, *sp_y = nullptr, *sp_z = nullptr, *sp_e = nullptr;

    tree->SetBranchAddress("px", &sim_px);
    tree->SetBranchAddress("py", &sim_py);
    tree->SetBranchAddress("pz", &sim_pz);
    tree->SetBranchAddress("purity", &sim_purity);
    tree->SetBranchAddress("completeness", &sim_completeness);
    tree->SetBranchAddress("x", &sp_x);
    tree->SetBranchAddress("y", &sp_y);
    tree->SetBranchAddress("z", &sp_z);
    tree->SetBranchAddress("e", &sp_e);

    TH1D* histx = new TH1D("histx", ";x (cm);Entries", 100, 0, 250);
    TH1D* histy = new TH1D("histy", ";y (cm);Entries", 100, -100, 100);
    TH1D* histz = new TH1D("histz", ";z (cm);Entries", 100, 0, 1000);

    TH1D* histpx = new TH1D("histpx", ";px (GeV/c);Entries", 30, 0, 1.5);
    TH1D* histpy = new TH1D("histpy", ";py (GeV/c);Entries", 30, 0, 1.5);
    TH1D* histpz = new TH1D("histpz", ";pz (GeV/c);Entries", 30, 0, 1.5);

    TH1D* histe = new TH1D("hitse", ";e (MeV/cm);Entries", 100, 0, 10); 

    TH1D* histpurity = new TH1D("histpurity", ";purity;Entries", 100, 0, 1);
    TH1D* histcompleteness = new TH1D("histcompleteness", ";completeness;Entries", 100, 0, 1);

    for (int it = 0; it < tree->GetEntries(); it++) {
        tree->GetEntry(it);

        for (double x : *sp_x) histx->Fill(x);
        for (double y : *sp_y) histy->Fill(y);
        for (double z : *sp_z) histz->Fill(z);
        for (double e : *sp_e) histe->Fill(e);

        histpx->Fill(sim_px);
        histpy->Fill(sim_py);
        histpz->Fill(sim_pz);

        std::cout << sim_purity << std::endl;
        std::cout << sim_completeness << std::endl;

        histpurity->Fill(sim_purity);
        histcompleteness->Fill(sim_completeness);
    }

    TCanvas* c1 = new TCanvas("c1", "", 800, 600);
    histx->Draw("HIST");
    c1->SaveAs("Plots/plot_x.png");
    delete c1;
    delete histx;

    TCanvas* c2 = new TCanvas("c2", "", 800, 600);
    histy->Draw("HIST");
    c2->SaveAs("Plots/plot_y.png");
    delete c2;
    delete histy;

    TCanvas* c3 = new TCanvas("c3", "", 800, 600);
    histz->Draw("HIST");
    c3->SaveAs("Plots/plot_z.png");
    delete c3;
    delete histz;

    TCanvas* c4 = new TCanvas("c4", "", 800, 600);
    histpx->Draw("HIST");
    c4->SaveAs("Plots/plot_px.png");
    delete c4;
    delete histpx;

    TCanvas* c5 = new TCanvas("c5", "", 800, 600);
    histpy->Draw("HIST");
    c5->SaveAs("Plots/plot_py.png");
    delete c5;
    delete histpy;

    TCanvas* c6 = new TCanvas("c6", "", 800, 600);
    histpz->Draw("HIST");
    c6->SaveAs("Plots/plot_pz.png");
    delete c6;
    delete histpz;

    TCanvas* c7 = new TCanvas("c7", "", 800, 600);
    histe->Draw("HIST");
    c7->SaveAs("Plots/plot_e.png");
    delete c7;
    delete histe;
    
    TCanvas* c8 = new TCanvas("c8", "", 800, 600);
    histpurity->Draw("HIST");
    c8->SaveAs("Plots/plot_purity.png");
    delete c8;
    delete histpurity;

    TCanvas* c9 = new TCanvas("c9", "", 800, 600);
    histcompleteness->Draw("HIST");
    c9->SaveAs("Plots/plot_completeness.png");
    delete c9;
    delete histcompleteness;

    file->Close();
    delete file;
}