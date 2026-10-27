// Calculation functions only. Configure and run with: root UCNDecayPlot.C
// Units: E,V,W in neV; time in s; S/volume in m^-1.
// Assumptions: isotropic, spatially uniform, one wall material,
// no gravity, no energy redistribution, Poisson collision rate v*S/(4*volume).
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <cmath>
#include <complex>
#include <vector>
#include <limits>
#include <algorithm>

namespace UCN {
// TMath has no neutron mass; obtain it from ROOT's particle database.
// PDG Mass() is in GeV/c^2; convert to kg with ROOT's SI constants.
double NeutronMassKg() {
  static const double mass = TDatabasePDG::Instance()->GetParticle(2112)->Mass()
      * 1e9 * TMath::Qe() / (TMath::C()*TMath::C());
  return mass;
}

// E in neV, speed in m/s. 1 neV = 1e-9 * Qe() joules.
double Speed(double E) {
  return std::sqrt(2*E*1e-9*TMath::Qe()/NeutronMassKg());
}

struct Quadrature {
  std::vector<double> x, w;
  explicit Quadrature(int n) : x(n), w(n) {
    const double pi = TMath::Pi();
    for (int i=0; i<(n+1)/2; ++i) {
      double z=std::cos(pi*(i+0.75)/(n+0.5)), pp=0;
      for (int k=0; k<100; ++k) {
        double p1=1, p2=0;
        for (int j=1; j<=n; ++j) {
          double p3=p2; p2=p1;
          p1=((2*j-1)*z*p2-(j-1)*p3)/j;
        }
        pp=n*(z*p1-p2)/(z*z-1);
        double dz=p1/pp; z-=dz;
        if (std::abs(dz)<1e-15) break;
      }
      x[i]=(1-z)/2; x[n-1-i]=(1+z)/2;
      w[i]=w[n-1-i]=1/((1-z*z)*pp*pp);
    }
  }
};

// Exact 1-|r|^2 for a semi-infinite wall U=V-iW, W>=0.
// q=sqrt(E_perp-V+iW); common sqrt(2m)/hbar cancels.
double LossPerBounce(double Eperp, double V, double W) {
  if (Eperp<=0) return 0;
  const double k=std::sqrt(Eperp);
  const auto q=std::sqrt(std::complex<double>(Eperp-V,W));
  // Algebraically 1-| (k-q)/(k+q) |^2; avoids cancellation.
  return 4*k*q.real()/std::norm(k+q);
}

// Collision flux weights angles as 2*cos(theta)*sin(theta)dtheta.
// Thus u=cos(theta)^2 is uniform and E_perp=E*u.
double MeanLoss(double E, double V, double W, int order=128) {
  const Quadrature q(order);
  double sum=0;
  for (int i=0;i<order;++i) sum+=q.w[i]*LossPerBounce(E*q.x[i],V,W);
  return sum;
}

// p: N0, V, W, Emax, S/volume, tau_beta, background, alpha (dN/dE proportional to E^alpha).
// tau_beta <= 0 disables beta decay. Emax must be <= V.
// constantMu < 0 selects optical-potential losses; otherwise use constantMu.
double DecayImpl(double *x, double *p, double constantMu) {
  const double t=x[0], V=p[1], W=p[2], emax=p[3], sv=p[4];
  if (t<0 || V<=0 || W<0 || emax<=0 || emax>V || sv<0 || p[7]<=-1 || constantMu>1)
    return std::numeric_limits<double>::quiet_NaN();
  static const Quadrature qe(128), qa(128);
  struct Cache {
    double V=-1,W=-1,emax=-1,sv=-1,alpha=-2,constantMu=-2;
    std::vector<double> rates, weights;
  };
  static Cache c; // For ordinary serial ROOT fits; not concurrent evaluation.
  if (V!=c.V || W!=c.W || emax!=c.emax || sv!=c.sv || p[7]!=c.alpha || constantMu!=c.constantMu) {
    c.V=V; c.W=W; c.emax=emax; c.sv=sv; c.alpha=p[7]; c.constantMu=constantMu;
    c.rates.resize(qe.x.size()); c.weights.resize(qe.x.size());
    double norm=0;
    for (unsigned i=0;i<qe.x.size();++i) {
      const double E=emax*qe.x[i];
      double mu=constantMu;
      if (constantMu<0) {
        mu=0;
        for (unsigned j=0;j<qa.x.size();++j)
          mu+=qa.w[j]*LossPerBounce(E*qa.x[j],V,W);
      }
      const double speed=Speed(E);
      c.rates[i]=speed*sv*mu/4;
      c.weights[i]=qe.w[i]*std::pow(qe.x[i],p[7]);
      norm+=c.weights[i];
    }
    for (auto &w:c.weights) w/=norm;
  }
  double survival=0;
  for (unsigned i=0;i<c.rates.size();++i)
    survival+=c.weights[i]*std::exp(-t*c.rates[i]);
  if (p[5]>0) survival*=std::exp(-t/p[5]);
  return p[0]*survival+p[6];
}
// Original 8-parameter interface, unchanged.
double Decay(double *x, double *p) { return DecayImpl(x,p,-1); }

// Same first 8 parameters as Decay(), plus p[8] = fixed loss per bounce.
// W does not enter the loss probability in this mode. Keep Emax <= V.
// The rate remains energy dependent: Gamma_wall(E)=mu0*Speed(E)*S/(4*volume).
double ConstantLossDecay(double *x, double *p) {
  if (!std::isfinite(p[8]) || p[8]<0 || p[8]>1)
    return std::numeric_limits<double>::quiet_NaN();
  return DecayImpl(x,p,p[8]);
}

// Angular mean to first order in W (weak absorption), retaining full E/V dependence.
// This is NOT a low-energy wall-loss model. Tiny-x series only avoids cancellation.
double MeanLossArcsin(double E,double V,double W) {
  if(E<0 || V<=0 || E>V || W<0)
    return std::numeric_limits<double>::quiet_NaN();
  const double u=E/V;
  if(u==0 || W==0) return 0;
  if(u<1e-4)
    return (4.0/3)*(W/V)*std::sqrt(u)*(1+u*(0.3+u*(9.0/56+u*5.0/48)));
  return 2*(W/V)*(std::asin(std::sqrt(u))/u-std::sqrt((1-u)/u));
}

// Fixed dN/dE ~ sqrt(E), y=v/vmax, normalized speed density 3*y^2.
// mu_bar(E) = c + MeanLossArcsin(E,V,W); additive small-loss model.
// Numerically integrate survival; no low-energy expansion in the fit model.
// p: N0, c_per_collision, W[neV], V, Emax, S/volume, tau_beta, B.
double MixedLossFit(double *x,double *p) {
  if(x[0]<0 || p[1]<0 || p[2]<0 || p[3]<=0 || p[4]<=0 || p[4]>p[3] || p[5]<0)
    return std::numeric_limits<double>::quiet_NaN();
  static const Quadrature quad(128);
  const double nu=Speed(p[4])*p[5]/4;
  double survival=0;
  for(unsigned i=0;i<quad.x.size();++i) {
    const double y=quad.x[i];
    const double mu=p[1]+MeanLossArcsin(p[4]*y*y,p[3],p[2]);
    survival+=3*quad.w[i]*y*y*std::exp(-x[0]*nu*y*mu);
  }
  return p[0]*(p[6]>0 ? std::exp(-x[0]/p[6]):1)*survival+p[7];
}

} // namespace UCN
