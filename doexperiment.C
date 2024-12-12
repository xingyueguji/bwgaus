#include "pesudoex.h"

void doexperiment()
{
    pesudoex* o1 = new pesudoex();

    o1->BWGaus = o1->dofftconvolution(o1->bw,o1->Gaus,o1->x);
    o1->GausGaus = o1->dofftconvolution(o1->bwGaus,o1->Gaus,o1->x);

    //cout << "The FWHM of BWGaus is" << o1->findFWHM(o1->BWGaus,o1->x) << endl;
    //cout << "The FWHM of GausGaus is" << o1->findFWHM(o1->GausGaus,o1->x) << endl;

    //cout << "The STD of BWGaus is" << o1->ComputeStdDevNumerical(o1->BWGaus,o1->x,10000) << endl;
    //cout << "The STD of GausGaus is" << o1->ComputeStdDevNumerical(o1->GausGaus,o1->x,10000) << endl;

    o1->plotPdf(o1->BWGaus,o1->bw,o1->x,"BWGaus");
    o1->plotPdf(o1->GausGaus,o1->bwGaus,o1->x,"GausGaus");


}