// Run from this directory: root UCNDecayPlot.C
#include "UCNDecay.C"
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TPad.h>
#include <iostream>
#include <TLatex.h>

// Speed spectrum for dN/dE proportional to sqrt(E).
// mu_bar(v) is a collision-flux angular average, NOT a convolution with dN/dv.
void DrawUCNVelocity(double N0, double V, double W, double Emax, double mu0, double S_over_volume,
                     const TString &outputDir) {
  const double vmax=UCN::Speed(Emax);
  auto canvas=new TCanvas("ucn_velocity_canvas", "UCN spectrum and angle-averaged wall loss",1200,1000);
  canvas->Divide(2,2);
  auto setupPad=[&](int i) {
    canvas->cd(i)->SetGrid();
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.14);
  };
  auto axes=[](TF1 *f) {
    f->GetXaxis()->SetTitleSize(0.042);
    f->GetYaxis()->SetTitleSize(0.04);
    f->GetYaxis()->SetTitleOffset(1.8);
    f->SetLineWidth(3);
    f->SetNpx(400);
    f->SetMinimum(0);
  };
  auto legend=[]() {
    auto l=new TLegend(0.23,0.75,0.92,0.88);
    l->SetBorderSize(0); l->SetTextSize(0.029); l->SetMargin(0.14);
    return l;
  };
  setupPad(1);
  auto spectrum=new TF1("ucn_speed_spectrum","3*[0]*x*x/([1]*[1]*[1])",0,vmax);
  spectrum->SetParNames("N0","vmax_m_per_s");
  spectrum->SetParameters(N0,vmax);
  spectrum->SetTitle("Initial speed distribution;Speed |v| [m/s];dN/d|v| [neutrons / (m/s)]");
  axes(spectrum);
  spectrum->SetLineColor(kBlue+1);
  spectrum->SetMaximum(1.3*3*N0/vmax);
  spectrum->Draw();
  auto l1=legend();
  l1->AddEntry(spectrum,"dN/dv #propto v^{2}","l");
  l1->AddEntry((TObject*)nullptr,Form("v_{max} = %.2f m/s",vmax),"");
  l1->Draw();

  setupPad(2);
  auto loss=new TF1("ucn_perpendicular_loss",[](double *x,double *p){
    const double Eperp=0.5*UCN::NeutronMassKg()*x[0]*x[0]/(1e-9*TMath::Qe());
    return UCN::LossPerBounce(Eperp,p[0],p[1]);
  },0,vmax,2);
  loss->SetParameters(V,W);
  loss->SetTitle("Loss per wall collision;Speed normal to the wall [m/s];Loss probability #mu(v_{n})");
  axes(loss);
  loss->SetLineColor(kRed+1);
  const double lossMax=1.3*std::max({loss->Eval(vmax),mu0,1e-8});
  loss->SetMaximum(lossMax);
  loss->Draw();
  auto constantNormal=new TF1("ucn_constant_normal_loss","[0]",0,vmax);
  constantNormal->SetParameter(0,mu0);
  constantNormal->SetLineColor(kGreen+2);
  constantNormal->SetLineStyle(2);
  constantNormal->SetLineWidth(3);
  constantNormal->Draw("SAME");
  auto l2=legend();
  l2->SetY1(0.70);
  l2->AddEntry(loss,"#mu(v_{n}) = 1 - R(v_{n})","l");
  l2->AddEntry((TObject*)nullptr,Form("V = %g, W = %g neV",V,W),"");
  l2->AddEntry(constantNormal,Form("Constant loss: #mu_{0} = %.0e",mu0),"l");
  l2->Draw();

  setupPad(3);
  auto mean=new TF1("ucn_speed_wall_loss",[](double *x,double *p){
    const double E=0.5*UCN::NeutronMassKg()*x[0]*x[0]/(1e-9*TMath::Qe());
    return UCN::MeanLoss(E,p[0],p[1]);
  },0,vmax,2);
  mean->SetParameters(V,W);
  mean->SetTitle("Angle-averaged wall loss;Speed |v| [m/s];Mean loss per bounce #bar{#mu}(v)");
  axes(mean);
  mean->SetLineColor(kBlack);
  mean->SetMaximum(lossMax); // Common vertical scale with the middle panel.
  mean->Draw();
  // Integral_0^v mu0 * 2*vn/v^2 dvn = mu0 (v>0).
  // At v=0 use the continuous limit; the wall collision frequency is zero.
  auto constantMean=new TF1("ucn_constant_mean_loss","[0]",0,vmax);
  constantMean->SetParameter(0,mu0);
  constantMean->SetLineColor(kGreen+2);
  constantMean->SetLineStyle(2);
  constantMean->SetLineWidth(3);
  constantMean->Draw("SAME");
  auto l3=legend();
  l3->SetY1(0.70);
  l3->AddEntry(mean,"#bar{#mu}(v): collision-flux average","l");
  l3->AddEntry((TObject*)nullptr,"Weight: 2v_{n}/v^{2},  0 #leq v_{n} #leq v","");
  l3->AddEntry(constantMean,Form("Averaged constant loss: %.0e",mu0),"l");
  l3->Draw();
  setupPad(4);
  // Initial wall-loss rate density: (dN/dv) * [v*S/(4*volume)] * mu_bar(v).
  // This is evaluated at t=0 and excludes beta decay.
  auto weighted=new TF1("ucn_weighted_speed",[=](double *x,double *){
    const double v=x[0];
    const double E=0.5*UCN::NeutronMassKg()*v*v/(1e-9*TMath::Qe());
    const double mu=UCN::MeanLoss(E,V,W);
    return (3*N0*v*v/(vmax*vmax*vmax))*(v*S_over_volume/4)*mu;
  },0,vmax,0);
  auto weightedConstant=new TF1("ucn_weighted_speed_constant",[=](double *x,double *){
    return (3*N0*x[0]*x[0]/(vmax*vmax*vmax))*(x[0]*S_over_volume/4)*mu0;
  },0,vmax,0);
  weighted->SetTitle("Initial wall-loss rate;Speed |v| [m/s];dN_{lost}/(dt dv) [neutrons / (s (m/s))]");
  axes(weighted); axes(weightedConstant);
  weighted->SetLineColor(kBlack);
  weightedConstant->SetLineColor(kGreen+2);
  weightedConstant->SetLineStyle(2);
  weighted->SetMaximum(1.35*std::max({weighted->Eval(vmax),weightedConstant->Eval(vmax),1e-8}));
  weighted->Draw(); weightedConstant->Draw("SAME");
  auto l4=legend();
  l4->SetY1(0.70);
  l4->AddEntry(weighted,"Velocity-dependent wall loss","l");
  l4->AddEntry(weightedConstant,Form("Constant loss: #mu_{0} = %.0e",mu0),"l");
  l4->AddEntry((TObject*)nullptr,"Includes wall collisions per second","");
  l4->Draw();
  std::cout << "Initial wall loss [neutrons/s]: optical=" << weighted->Integral(0,vmax)
            << ", constant=" << weightedConstant->Integral(0,vmax) << std::endl;
  canvas->Modified(); canvas->Update();
  canvas->SaveAs(outputDir+"/UCNVelocityPlot.png");
  canvas->SaveAs(outputDir+"/UCNVelocityPlot.pdf");
  std::cout<<"Speed spectrum: alpha=0.5, vmax="<<vmax<<" m/s, integral="
           <<spectrum->Integral(0,vmax)<<", eta="<<W/V<<std::endl;
}


void DrawUCNMixedScans(const double *base, const std::vector<double>& muScan,
 const std::vector<double>& wScan,double tmin,double tmax,double ymin,double ymax,
 int npx,const TString &outputDir) {
  auto make=[&](const char *name) {
    auto f=new TF1(name,UCN::MixedLossFit,tmin,tmax,8);
    f->SetParNames("N0","c_per_collision","W_neV","V_neV","Emax_neV",
                   "S_over_volume_mInv","tau_beta_s","background");
    f->SetParameters(base);f->SetNpx(npx);
    f->SetParLimits(0,0,1e12);f->SetParLimits(1,0,0.01);
    f->SetParLimits(2,0,0.01*base[3]);
    for(int j=3;j<8;++j) f->FixParameter(j,base[j]);
    return f;
  };
  make("fit_mixed_loss");

  const std::vector<double> scans[]={muScan,wScan};
  const int colors[]={kBlue+1,kRed+1,kMagenta+1,kOrange+7};
  for(int panel=0;panel<2;++panel) {
    auto c=new TCanvas(Form("ucn_mixed_fit_canvas_%d",panel),
      panel==0 ? "Constant loss component" : "Optical loss component",1000,700);
    c->SetLogy();c->SetGrid();c->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    auto l=new TLegend(0.18,0.18,0.64,0.36);l->SetBorderSize(0);l->SetTextSize(0.031);
    for(unsigned j=0;j<scans[panel].size();++j) {
      auto f=make(Form("mixed_scan_%d_%u",panel,j));
      f->SetParameter(panel+1,scans[panel][j]);
      f->SetTitle(Form("%s;Storage time [s];Surviving UCN",panel==0 ? "Constant loss component":"Optical loss component"));
      f->SetMinimum(ymin);f->SetMaximum(ymax);f->GetYaxis()->SetTitleOffset(1.5);
      const bool baseline = scans[panel][j] == base[panel+1];
      f->SetLineColor(baseline ? kBlack : colors[j%4]);
      f->SetLineWidth(3);f->SetLineStyle(baseline ? 1 : 2+j);
      f->Draw(j?"SAME":"");
      l->AddEntry(f,panel==0?Form("c = %.1e",scans[panel][j]):Form("W = %.3g neV",scans[panel][j]),"l");
    }
    l->Draw();auto note=new TLatex();note->SetNDC();note->SetTextSize(0.029);
    note->DrawLatex(0.18,0.49,"Fixed spectrum: dN/dE #propto E^{1/2}");
    note->DrawLatex(0.18,0.44,panel==0?Form("Fixed W = %.3g neV",base[2]):Form("Fixed c = %.1e",base[1]));
    const TString stem=panel==0 ? "/UCNMixedMuScan" : "/UCNMixedWScan";
    c->Modified();c->Update();
    c->SaveAs(outputDir+stem+".png");c->SaveAs(outputDir+stem+".pdf");
  }
  std::cout<<"Ready: graph->Fit(\"fit_mixed_loss\",\"RS\"); free N0, c_per_collision, W; spectrum fixed.\n";
}

void UCNDecayPlot() {
  // All model inputs are specified here. Values are illustrative.
  const double N0 = 10000;
  const double V = 200.0;           // wall potential [neV]
  const double W = 0.02;            // imaginary potential magnitude [neV]
  const double Emax = 150.0;        // spectral cutoff [neV], 0 < Emax <= V
  const double S_over_volume = 30; // inner surface / internal volume [m^-1]
  const double tau_beta = 880.0;   // beta lifetime [s]; <= 0 disables beta decay
  const double B = 0.0;            // no constant residual UCN population

  // mu0 is a probability, NOT the optical ratio W/V.
  const double mu0 = 1.0e-4;       // fixed loss per bounce, illustrative value

  // Display options: true = show, false = hide (including legend entry).
  const bool showVelocityPlots = true; // speed spectrum and wall-loss comparison
  const bool showExponential = true;
  const bool showPractical = true;
  const bool showConstantLoss = true;
  const int visibleCount = showExponential + showPractical + showConstantLoss;
  if (showVelocityPlots) DrawUCNVelocity(N0,V,W,Emax,mu0,S_over_volume,gSystem->DirName(__FILE__));
  if (visibleCount == 0) {
    std::cout << "No curves selected; enable at least one display option. No decay files saved." << std::endl;
    return;
  }

  const double alpha = 0.5;        // dN/dE proportional to sqrt(E); N(<E) proportional to E^(3/2)

  // Plot range and numerical initial-slope step.
  const double tmin = 0.0, tmax = 2000.0; // [s]
  const int npx = 500;
  const bool showMixedFitScans = true;
  const double mixedMu = 0.0; // baseline matches the black reference curve
  const double mixedW = W; // W [neV] for the arcsin model
  const std::vector<double> scanMu = {0.0, 0.5e-4, 1.0e-4};
  const std::vector<double> scanW = {0.01, 0.02, 0.04};

  const double h = 0.01; // [s]

  auto f = new TF1("f", UCN::Decay, tmin, tmax, 8);
  f->SetParNames("N0", "V_neV", "W_neV", "Emax_neV",
                 "S_over_volume_mInv", "tau_beta_s", "background", "alpha");
  f->SetParameters(N0, V, W, Emax, S_over_volume, tau_beta, B, alpha);
  // Fitting setup: spectrum fixed; N0 and W free.
  f->SetParLimits(7, 0.0, 5.0);
  f->SetParLimits(0, 0, 1e12);
  f->SetParLimits(2, 0, 2);
  for (int i : {1, 3, 4, 5, 6, 7})
    f->FixParameter(i, f->GetParameter(i));
  f->SetNpx(npx);
  f->SetTitle("UCN storage;Storage time [s];Surviving UCN");

  // Match initial population and initial logarithmic decay rate.
  // Second-order forward difference avoids evaluation at negative time.
  const double background = f->GetParameter(6);
  const double n0 = f->Eval(0) - background;
  const double log1 = std::log((f->Eval(h) - background) / n0);
  const double log2 = std::log((f->Eval(2*h) - background) / n0);
  const double rate = -(4*log1 - log2)/(2*h);
  auto exponential = new TF1("single_exp", "[0]*exp(-[1]*x)+[2]",
                             f->GetXmin(), f->GetXmax());
  exponential->SetParNames("N0", "decay_rate_sInv", "background");
  exponential->SetParameters(n0, rate, background);
  exponential->SetNpx(npx);

  // Reference decay curves.
  auto canvas = new TCanvas("ucn_canvas", "UCN decay comparison", 1000, 700);
  canvas->SetGrid();
  canvas->SetLogy();
  canvas->SetLeftMargin(0.15);
  canvas->SetBottomMargin(0.14);
  auto legend = new TLegend(0.15, 0.15, 0.61, 0.15+0.0575*visibleCount);
  legend->SetBorderSize(0);
  f->SetLineColor(kBlack);
  f->SetLineWidth(3);
  exponential->SetLineColor(kRed+1);
  exponential->SetLineStyle(2);
  exponential->SetLineWidth(3);
  std::vector<TF1*> curves;
  if (showExponential) {
    legend->AddEntry(exponential, "Single exponential", "l");
    curves.push_back(exponential);
  }
  if (showPractical) {
    legend->AddEntry(f, "Numerical wall-loss model", "l");
    curves.push_back(f);
  }
  if (showConstantLoss) {
    auto constantLoss=new TF1("constant_loss",UCN::ConstantLossDecay,tmin,tmax,9);
    for (int j=0;j<8;++j) {
      constantLoss->SetParName(j,f->GetParName(j));
      constantLoss->SetParameter(j,f->GetParameter(j));
    }
    constantLoss->SetParName(8,"mu0_per_bounce");
    constantLoss->FixParameter(8,mu0);
    // W is inactive here. Keep geometry, spectrum and beta lifetime fixed.
    for (int j : {1,2,3,4,5,6,7}) constantLoss->FixParameter(j,constantLoss->GetParameter(j));
    constantLoss->SetParLimits(0,0,1e12);
    constantLoss->SetNpx(npx);
    constantLoss->SetLineColor(kGreen+2);
    constantLoss->SetLineStyle(9);
    constantLoss->SetLineWidth(3);
    legend->AddEntry(constantLoss,Form("Constant loss per bounce (#mu_{0} = %.1e)",mu0),"l");
    curves.push_back(constantLoss);
  }
  double ymin=curves.front()->Eval(tmax);
  for (auto curve : curves) ymin=std::min(ymin,curve->Eval(tmax));
  curves.front()->SetTitle("UCN storage;Storage time [s];Surviving UCN");
  curves.front()->GetYaxis()->SetTitleOffset(1.5);
  curves.front()->SetMinimum(std::max(0.5*ymin,1e-12));
  curves.front()->SetMaximum(2*N0);
  for (unsigned i=0;i<curves.size();++i) curves[i]->Draw(i ? "SAME" : "");
  legend->Draw();
  canvas->Modified();
  canvas->Update();

  // Save next to this macro, also when invoked from another directory.
  const TString outputDir=gSystem->DirName(__FILE__);
  canvas->SaveAs(outputDir+"/UCNDecayPlot.png");
  canvas->SaveAs(outputDir+"/UCNDecayPlot.pdf");
  if (showMixedFitScans) {
    const double base[]={N0,mixedMu,mixedW,V,Emax,S_over_volume,tau_beta,B};
    DrawUCNMixedScans(base,scanMu,scanW,tmin,tmax,
                        std::max(0.5*ymin,1e-12),2*N0,npx,outputDir);
  }
  if (showExponential && rate>0)
    std::cout << "Single exponential: same N(0) and initial slope; tau = "
              << 1/rate << " s" << std::endl;
}
