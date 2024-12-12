#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <TCanvas.h>
#include <TStyle.h>

class pesudoex
{

    // 1. first do a pdf + pdf convolution, see the difference in width
    // 2. second do a smearing simulating the template, then fit the histogram.
    // 3. third do a varible binning, with limited statistics.

public:
    pesudoex();
    ~pesudoex();

    RooFFTConvPdf *dofftconvolution(RooAbsPdf *pdf1, RooAbsPdf *pdf2, RooRealVar *x);
    double findMaxValue(RooAbsPdf *pdf, RooRealVar *x);
    double findFWHM(RooAbsPdf *pdf, RooRealVar *x);
    void plotPdf(RooAbsPdf *pdf, RooAbsPdf *pdf2, RooRealVar *x, const char *canvasName = "pdfCanvas");
    double ComputeStdDevNumerical(RooAbsPdf *pdf, RooRealVar *var, int nSteps = 1000);
    double ComputeMeanNumerical(RooAbsPdf *pdf, RooRealVar *var, int nSteps = 1000);

    // Member functions

    // Data members

    RooRealVar *x;
    // This is Breit_Wigner and Gaus like Breit_wigner
    RooRealVar *bwmean;
    RooRealVar *width;
    RooBreitWigner *bw;

    RooRealVar *bwgausmean;
    RooRealVar *bwgaussig;
    RooGaussian *bwGaus;

    // This is Gaussian
    RooRealVar *Gausmean;
    RooRealVar *Gaussigma;
    RooGaussian *Gaus;

    RooFFTConvPdf *GausGaus;
    RooFFTConvPdf *BWGaus;
};

pesudoex::pesudoex()
{
    x = new RooRealVar("roomass", "roomass", 60, 120);

    x->setBinning(RooBinning(10000, 60, 120), "cache");
    bwmean = new RooRealVar("bwmean", "bwmean", 91);
    width = new RooRealVar("width", "width", 2.5);
    bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width);

    bwgausmean = new RooRealVar("bwgausmean", "bwgausmean", 91);
    bwgaussig = new RooRealVar("bwgaussig", "bwgaussig", 2.5 / 2.355);
    bwGaus = new RooGaussian("bwGaus", "bwGaus", *x, *bwgausmean, *bwgaussig);

    Gausmean = new RooRealVar("Gausmean", "Gausmean", 0);
    Gaussigma = new RooRealVar("Gaussigma", "Gaussigma", 0.6);
    Gaus = new RooGaussian("Gaus", "Gaus", *x, *Gausmean, *Gaussigma);
}

double pesudoex::ComputeMeanNumerical(RooAbsPdf *pdf, RooRealVar *var, int nSteps = 1000)
{
    if (!pdf || !var || nSteps <= 0)
    {
        std::cerr << "Invalid input: PDF, variable, or number of steps." << std::endl;
        return -1.0;
    }

    // Get the range of the variable
    double xMin = var->getMin();
    double xMax = var->getMax();
    double stepSize = (xMax - xMin) / nSteps;

    // Initialize mean calculation
    double mean = 0.0;
    double norm = 0.0;

    // Iterate over the range of the variable
    for (int i = 0; i <= nSteps; ++i)
    {
        double x = xMin + i * stepSize;
        var->setVal(x);                    // Set the value of the variable
        double pdfVal = pdf->getVal(*var); // Evaluate the PDF at this value
        double weight = pdfVal * stepSize; // Weight for this step
        mean += weight * x;                // Accumulate weighted value of x
        norm += weight;                    // Accumulate normalization
    }

    // Normalize mean
    mean /= norm;

    return mean;
}
double pesudoex::ComputeStdDevNumerical(RooAbsPdf *pdf, RooRealVar *var, int nSteps = 1000)
{
    if (!pdf || !var || nSteps <= 0)
    {
        std::cerr << "Invalid input: PDF, variable, or number of steps." << std::endl;
        return -1.0;
    }

    // Get the range of the variable
    double xMin = var->getMin();
    double xMax = var->getMax();
    double stepSize = (xMax - xMin) / nSteps;

    // Initialize variance calculation
    double variance = 0.0;
    double norm = 0.0;

    double mean = ComputeMeanNumerical(pdf,var,nSteps);

    // Iterate over the range of the variable
    for (int i = 0; i <= nSteps; ++i)
    {
        double x = xMin + i * stepSize;
        var->setVal(x);                               // Set the value of the variable
        double pdfVal = pdf->getVal(*var);            // Evaluate the PDF at this value
        double weight = pdfVal * stepSize;            // Weight for this step
        variance += weight * (x - mean) * (x - mean); // Accumulate variance
        norm += weight;                               // Accumulate normalization
    }

    // Normalize variance
    variance /= norm;

    // Standard deviation is the square root of variance
    return std::sqrt(variance);
}

RooFFTConvPdf *pesudoex::dofftconvolution(RooAbsPdf *pdf1, RooAbsPdf *pdf2, RooRealVar *x)
{
    RooFFTConvPdf *temp = new RooFFTConvPdf("temp", "temp", *x, *pdf1, *pdf2);

    return temp;
}

double pesudoex::findMaxValue(RooAbsPdf *pdf, RooRealVar *x)
{
    double maxVal = -1.0;      // Initialize with a small value
    double xMax = x->getMin(); // Initialize max position

    double step = (x->getMax() - x->getMin()) / 1000000.0; // Adjust resolution
    for (double xx = x->getMin(); xx <= x->getMax(); xx += step)
    {
        x->setVal(xx);              // Set the value of x
        double val = pdf->getVal(); // Evaluate the PDF
        if (val > maxVal)
        {
            maxVal = val; // Update max value
            xMax = xx;    // Update max position
        }
    }

    return maxVal;
}

double pesudoex::findFWHM(RooAbsPdf *pdf, RooRealVar *x)
{
    // Step 1: Find the maximum value
    double maxVal = this->findMaxValue(pdf, x);

    // Step 2: Define half maximum
    double halfMax = maxVal / 2.0;

    // Step 3: Scan for FWHM boundaries
    double xLow = x->getMin();
    double xHigh = x->getMax();
    double step = (xHigh - xLow) / 1000000.0; // Adjust resolution as needed

    bool foundLow = false, foundHigh = false;
    double lowerBound = xLow, upperBound = xHigh;

    for (double xx = xLow; xx <= xHigh; xx += step)
    {
        x->setVal(xx);
        double val = pdf->getVal();
        if (!foundLow && val >= halfMax)
        {
            lowerBound = xx;
            foundLow = true;
        }
        if (foundLow && val < halfMax)
        {
            upperBound = xx;
            foundHigh = true;
            break;
        }
    }

    // Step 4: Calculate FWHM
    return upperBound - lowerBound;
}

void pesudoex::plotPdf(RooAbsPdf *pdf, RooAbsPdf *pdf2, RooRealVar *x, const char *canvasName = "pdfCanvas")
{
    // Create a canvas
    TCanvas *canvas = new TCanvas(canvasName, "PDF Plot", 800, 600);

    // Set up general canvas cosmetics
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);

    // Create a RooPlot frame for the observable
    RooPlot *frame = x->frame();
    frame->SetTitle("PDF Plot");
    frame->GetXaxis()->SetTitle("x");      // X-axis title
    frame->GetXaxis()->SetTitleSize(0.05); // X-axis title size
    frame->GetXaxis()->SetTitleFont(42);   // X-axis title font (Helvetica)
    frame->GetXaxis()->SetLabelSize(0.04); // X-axis label size
    frame->GetXaxis()->SetLabelFont(42);   // X-axis label font (Helvetica)

    frame->GetYaxis()->SetTitle("Probability Density");
    frame->GetYaxis()->SetTitleSize(0.05); // Y-axis title size
    frame->GetYaxis()->SetTitleFont(42);   // Y-axis title font (Helvetica)
    frame->GetYaxis()->SetLabelSize(0.04); // Y-axis label size
    frame->GetYaxis()->SetLabelFont(42);   // Y-axis label font (Helvetica)

    // Plot the PDF on the frame
    pdf->plotOn(frame, RooFit::LineColor(kBlue), RooFit::LineWidth(2));
    pdf2->plotOn(frame, RooFit::LineColor(kRed), RooFit::LineWidth(2));

    // Draw the frame
    frame->Draw();

    TPaveText *pt = new TPaveText(0.2, 0.65, 0.35, 0.85, "brNDC"); // normalized coordinates
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextAlign(12); // Align left and vertically centered
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    pt->SetMargin(0.02);
    double FWHM = this->findFWHM(pdf, x);
    double FWHM2 = this->findFWHM(pdf2, x);
    double STD = this->ComputeStdDevNumerical(pdf,x,1000000);
    double STD2 = this->ComputeStdDevNumerical(pdf2,x,1000000);
    pt->AddText(Form("FWHM = %.4f", FWHM));
    pt->AddText(Form("STD = %.4f", STD));
    pt->AddText(Form("Assumption STD + STD = %.4f", sqrt(pow(STD2,2) + pow(0.6,2)) ));
    pt->AddText(Form("Assumption STD (Gamma) + STD = %.4f", sqrt(pow(FWHM2,2) + pow(0.6,2)) ));
    pt->AddText(Form("Assumption FWHM from STD (Gamma) + STD = %.4f", 2.355* sqrt(pow(FWHM2,2) + pow(0.6,2)) ));

    pt->Draw();

    // Optional: Draw cosmetics
    gStyle->SetOptStat(0); // Disable statistics box
    gPad->Modified();
    gPad->Update();

    canvas->SaveAs(Form("%s.png", canvasName));
}