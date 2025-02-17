#include "createratio.C"

void run_ratio()
{
  // format: Name, shiftlow, shifthigh, smearlow, smearhigh
  /*createratio("eta_0_10.root", -0.4, -0.08, -0.25, 0.3);
  createratio("eta_10_20.root", -0.25, 0.15, -0.05, 0.65);
  createratio("eta_20_30.root", -0.4, 0.075, -0.15, 0.65);
  createratio("eta_30_100.root", -0.45, -0.05, -0.3, 0.5);
  createratio("eta_0_100.root", -0.28, -0.07, 0.03, 0.36);

  createratio("raw_0_10.root", -0.24, -0.02, 0.05, 0.4);
  createratio("raw_10_20.root", -0.32, -0.06, 0.175, 0.6);
  createratio("raw_20_30.root", -0.32, -0.04, -0.1, 0.4);
  createratio("raw_30_100.root", -0.24, 0.02, 0.2, 0.7);
  createratio("raw_0_100.root", -0.22, -0.08, 0.22, 0.45);*/

  createratio("pp_eta_bk.root", -0.14, -0.06, 0.2, 0.32);
  createratio("pp_raw_bk.root", -0.105, -0.05, 0.2, 0.3);
  createratio("pp_eta_nobk.root", -0.15, -0.06, 0.23, 0.37);
  createratio("pp_raw_nobk.root", -0.09, -0.05, 0.28, 0.36);

  createratio("pp_eta_mass_range.root", -0.14, -0.07, 0.13, 0.25);
  createratio("pp_raw_mass_range.root", -0.1, -0.05, 0.16, 0.22);

}