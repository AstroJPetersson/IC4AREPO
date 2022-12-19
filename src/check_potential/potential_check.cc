#include "./common.h"
#include "./potential.h"

int main(int argc, char **argv) {
        
    /*Internal units*/

    double InternalMass_in_cgs     = 1.988e33;
    double InternalVelocity_in_cgs = 1.e5;
    double InternalLength_in_cgs   = 3.086e20;
    
    
    double G = 6.67e-8 * InternalMass_in_cgs / InternalVelocity_in_cgs / InternalVelocity_in_cgs / InternalLength_in_cgs;

    ExpBarPotential* Barpot;
    ExpDiskPotential* Diskpot1;
    ExpDiskPotential* Diskpot2;
    DPLPotential* Halopot;
    ModifiedMcMillanBulgePotential* Bulgepot;

    /*Initializing bar potential*/
    double rminbar = 0., rmaxbar = 0.;
    int nrbar = 500, npolybar = 6, ngaussbar = 8;
    double rho0bar = 5.e6, a0bar = 7.5, qbar = 0.5;

    ExpBarDensity   rhobar_axi(rho0bar*qbar*qbar, a0bar, 1., rminbar, rmaxbar);
    Barpot = new ExpBarPotential(G, rhobar_axi, nrbar, npolybar, ngaussbar);

    /*Initializing Bulge potential*/
    double rminbulge = 0., rmaxbulge = 0.;
    int nrbulge = 500, npolybulge = 6, ngaussbulge = 8;
    double rho0bulge = 8.0e5, a0bulge=10.0, acutbulge = 10.0;
    double alfbulge = 1.7, qbulge = 0.5;

    ModifiedMcMillanBulgeDensity   rhobulge(rho0bulge, alfbulge, a0bulge, acutbulge, qbulge, rminbulge, rmaxbulge);
    Bulgepot = new ModifiedMcMillanBulgePotential(G, rhobulge, nrbulge, npolybulge, ngaussbulge);

    /*Initializing exponential disk potential
     *Thick stellar disk*/
    double rmind1 = 0., rmaxd1 = 0.;
    int nrd1 = 500, npolyd1 = 20, ngaussd1 = 8;
    double Sigma0d1 = 1.74e6, Rd1 = 30.2, zd1 = 9.;

    ExpDiskDensity   rhodisk1(Sigma0d1, zd1, Rd1, rmind1, rmaxd1);
    Diskpot1 = new ExpDiskPotential(G, rhodisk1, nrd1, npolyd1, ngaussd1);

    /*Thin stellar disk*/
    double rmind2 = 0., rmaxd2 = 0.;
    int nrd2 = 500, npolyd2 = 20, ngaussd2 = 8;
    double Sigma0d2 = 8.5e6, Rd2 = 25.0, zd2 = 3.0;

    ExpDiskDensity   rhodisk2(Sigma0d2, zd2, Rd2, rmind2, rmaxd2);
    Diskpot2 = new ExpDiskPotential(G, rhodisk2, nrd2, npolyd2, ngaussd2);

    /*Initializing Halo potential*/
    double rminh = 0., rmaxh = 0.;
    int nrh = 500, npolyh = 6, ngaussh = 8;
    double rho0h = 8.11e3, ah = 196., alphah = 1.0, betah = 3.0;

    DPLDensity   rhohalo(rho0h, ah, alphah, betah, rminh, rmaxh);
    Halopot = new DPLPotential(G, rhohalo, nrh, npolyh, ngaussh);

    
    
   
    double Dx = 1000.;//internal units
    double dx = 0.5; //internal units 

    int nx = Dx/dx;
      
    double x0 = 0.;
    
    stringstream outfile_vc;
    outfile_vc << "./vc.txt";
    ofstream fout(outfile_vc.str().c_str());
    
    fout << x0 << "\t" << 0. << endl;
      
    for(int ix=1; ix<nx; ix++)
      {
        double x = x0 + ix*dx;
        Vec3 X(x, 0, 0);
        Vec3 dPhidx = Barpot->dPhidx(X) + Bulgepot->dPhidx(X) + Diskpot1->dPhidx(X) + Diskpot2->dPhidx(X) + Halopot->dPhidx(X);
        double Vc = sqrt(dPhidx[0] * x);
        fout << x << "\t" << Vc << endl;
      }
    
    cout << "Vc saved to vc.txt" << endl << endl;
}
