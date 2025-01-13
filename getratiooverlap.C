#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>

void getratiooverlap(TString data = "eta_0_10.root", TString theory = "0-10_theory.root")
{
    const int nbins_cent = 11;
    const int nbins_smear = 21;
    const int nbins_mass_shift = 21;

    double eta_mass_shift_array_low[nbins_cent] = {-0.45, -0.3, -0.45, -0.5, 0, 0, 0, 0, 0, 0, -0.35};
    double eta_mass_shift_array_high[nbins_cent] = {-0.05, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0.0};
    double eta_mass_smear_array_low[nbins_cent] = {-0.25, -0.05, -0.15, -0.3, 0, 0, 0, 0, 0, 0, 0};
    double eta_mass_smear_array_high[nbins_cent] = {0.3, 0.65, 0.65, 0.5, 0, 0, 0, 0, 0, 0, 0.35};

    double raw_mass_shift_array_low[nbins_cent] = {-0.3, -0.35, -0.4, -0.3, 0, 0, 0, 0, 0, 0, -0.25};
    double raw_mass_shift_array_high[nbins_cent] = {0.05, 0.0, 0.0, 0.1, 0, 0, 0, 0, 0, 0, 0.0};
    double raw_mass_smear_array_low[nbins_cent] = {0.05, 0.15, -0.1, 0.2, 0, 0, 0, 0, 0, 0, 0.2};
    double raw_mass_smear_array_high[nbins_cent] = {0.4, 0.6, 0.4, 0.7, 0, 0, 0, 0, 0, 0, 0.45};

    TFile *f_theory = new TFile(theory, "READ");
    TFile *f_pesudoexp = new TFile(data, "READ");

    TF1 *tf_theory = (TF1 *)f_theory->Get("Compare_shift_0_smear_0");
    tf_theory->SetLineColor(kRed);
    tf_theory->SetLineWidth(2);

    TF1 *tf_pesudoexp[nbins_mass_shift][nbins_smear];

    TCanvas *c1 = new TCanvas("", "", 800, 600);
    c1->cd();

    for (int i = 0; i < nbins_mass_shift; i++)
    {
        for (int j = 0; j < nbins_smear; j++)
        {
            tf_pesudoexp[i][j] = (TF1 *)f_pesudoexp->Get(Form("Compare_shift_%i_smear_%i", i, j));
            tf_pesudoexp[i][j]->SetLineColor(kBlue);
            tf_pesudoexp[i][j]->SetLineColorAlpha(kBlue, 0.05);
        }
    }
    // tf_theory->GetYaxis()->SetRangeUser(0.6,1.6);
    tf_theory->SetTitle("");
    tf_theory->GetYaxis()->SetTitle("Ratio");
    tf_theory->GetXaxis()->SetTitle("Mass (GeV)");

    tf_theory->Draw();

    for (int i = 0; i < nbins_mass_shift; i++)
    {
        for (int j = 0; j < nbins_smear; j++)
        {
            tf_pesudoexp[i][j]->Draw("SAME");
        }
    }

    tf_theory->Draw("SAME");

    auto legend = new TLegend(0.2, 0.2, 0.4, 0.4);
    legend->AddEntry(tf_theory, "Theory (Red)", "l");
    legend->AddEntry(tf_pesudoexp[0][0], "Template (Shaded)", "l");
    legend->Draw();


    TString newdata = data;
    TString newtheory = theory;

    newdata.Remove(newdata.Last('.'));
    newtheory.Remove(newtheory.Last('.'));
    
    TString ratioplotname = newdata + "_vs_" + newtheory;

    c1->SaveAs(Form("./overlapplot/%s.png", ratioplotname.Data()));
}