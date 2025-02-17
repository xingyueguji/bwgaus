#include "pesudoex.h"

void createratio(TString filename = "20-30_theory.root", double shiftlowbin = -0.2451, double shifthighbin = -0.2451, double smearlowbin = 0.3062, double smearhighbin = 0.3062)
{

    /*ZInvMassTheory = {-0.06287, -0.1768, -0.2451, -0.3271, -0.3362, -0.3294, -0.2702};
    ZWidthTheory = {0.08064, 0.2219, 0.3062, 0.4018, 0.4246, 0.4155, 0.354};

    Each of these numbers in the arrays above are organized in the following centrality scheme : 
    (0 - 10), (10 - 20), (20 - 30), (30 - 40), (40 - 50), (50 - 60) and (60 - 70).*/

    const int nbins_mass_shift = 42;
    const int nbins_smear = 42;
    const int nbins_cent = 11;

    double shift_bin[nbins_mass_shift];
    double smear_bin[nbins_smear];

    double shiftamount = (shifthighbin - shiftlowbin) / (nbins_mass_shift - 1); // 21 ticks from -0.5 to 0.5, 20 divisions
    double smearedamount = (smearhighbin - smearlowbin) / (nbins_smear - 1);    // 21 ticks from 0 to 0.2, 20 divisions.

    cout << "Shift amount is " << shiftamount << " Smear amount is " << smearedamount << endl;

    for (int i = 0; i < nbins_mass_shift; i++)
    {
        shift_bin[i] = shiftlowbin + i * shiftamount;
        if (fabs(shift_bin[i]) < 1e-10)
            shift_bin[i] = 0.0;
        cout << "shift bin " << i << " is " << shift_bin[i] << endl;
    }
    for (int j = 0; j < nbins_smear; j++)
    {
        smear_bin[j] = smearlowbin + j * smearedamount;
        if (fabs(smear_bin[j]) < 1e-10)
            smear_bin[j] = 0.0;
        cout << "smear bin " << j << " is " << smear_bin[j] << endl;
    }

    TNamed *namedString = new TNamed("Description", Form("This is shift %.5f_%.5f smear %.5f_%.5f", shiftlowbin, shifthighbin, smearlowbin, smearhighbin));

    for (int i = 0; i < nbins_mass_shift; i++)
    {
        for (int j = 0; j < nbins_smear; j++)
        {
            pesudoex *o1 = new pesudoex(91.1876, 2.4955, 91.1876 + shift_bin[i], 2.4955 + smear_bin[j]);
            o1->plotPdf(o1->bw, o1->bw2, o1->x, namedString, Form("Compare_shift_%i_smear_%i", i, j), filename);
        }
    }
}