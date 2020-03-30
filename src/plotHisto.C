// ROOT macro file for plotting example B4 histograms 
// 
// Can be run from ROOT session:
// root[0] .x plotHisto.C

{

//................
// Define Style
//
TStyle *nnbarStyle = new TStyle("NNbar","NNbar style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  nnbarStyle->SetFrameBorderMode(icol);
  nnbarStyle->SetFrameFillColor(icol);
  nnbarStyle->SetCanvasBorderMode(icol);
  nnbarStyle->SetCanvasColor(icol);
  nnbarStyle->SetPadBorderMode(icol);
  nnbarStyle->SetPadColor(icol);
  nnbarStyle->SetStatColor(icol);
  //nnbarStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  //Legend
  Int_t ileg=0;
  nnbarStyle->SetLegendBorderSize(ileg);
  nnbarStyle->SetLegendFillColor(ileg);
  nnbarStyle->SetLegendTextSize(0.045);
  nnbarStyle->SetLegendFont(42);

  // set the paper & margin sizes
  nnbarStyle->SetPaperSize(20,26);

  // set margin sizes
  nnbarStyle->SetPadTopMargin(0.07);
  nnbarStyle->SetPadRightMargin(0.05);
  nnbarStyle->SetPadBottomMargin(0.16);
  nnbarStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  nnbarStyle->SetTitleXOffset(1.5);
  nnbarStyle->SetTitleYOffset(1.5);


  gROOT->Reset();
  gROOT->SetStyle("NNbar");
  
//......................................................

  // Draw histos filled by Geant4 simulation 
  //  
  //
  //  Eabs
  //  Egap
  //  Labs
  //  NumCerenkov
  //  PhotonTime
  //  DecayTime
  //  PhotonPos 

  // Open file filled by Geant4 simulation 
  TFile f("B4.root");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  c1->Divide(2,2);
  
  // Draw Eabs histogram in the pad 1
  c1->cd(1);
  TH1D* hist1 = (TH1D*)f.Get("Eabs");
  hist1->Draw("HIST");
  
  // Draw Labs histogram in the pad 2
  c1->cd(2);
  TH1D* hist2 = (TH1D*)f.Get("Egap");
  hist2->Draw("HIST");
  
  // Draw Egap histogram in the pad 3
  // with logaritmic scale for y
  TH1D* hist3 = (TH1D*)f.Get("NumCerenkov");
  c1->cd(3);
  gPad->SetLogy(1);
  hist3->Draw("HIST");
  
  // Draw Lgap histogram in the pad 4
  // with logaritmic scale for y
  c1->cd(4);
  gPad->SetLogy(1);
  TH1D* hist4 = (TH1D*)f.Get("PhotonTime");
  hist4->Draw("HIST");

  // set label offset
  nnbarStyle->SetLabelOffset(0.01,"xyz");

  // use large fonts
  Int_t font=42; //  helvetica-medium-r-normal "Arial"
  Double_t tsize=0.055;
  nnbarStyle->SetTextFont(font);

  nnbarStyle->SetTextSize(tsize);
  nnbarStyle->SetLabelFont(font,"x");
  nnbarStyle->SetTitleFont(font,"x");
  nnbarStyle->SetLabelFont(font,"y");
  nnbarStyle->SetTitleFont(font,"y");
  nnbarStyle->SetLabelFont(font,"z");
  nnbarStyle->SetTitleFont(font,"z");

  nnbarStyle->SetLabelSize(tsize,"x");
  nnbarStyle->SetTitleSize(tsize,"x");
  nnbarStyle->SetLabelSize(tsize,"y");
  nnbarStyle->SetTitleSize(tsize,"y");
  nnbarStyle->SetLabelSize(tsize,"z");
  nnbarStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  nnbarStyle->SetMarkerStyle(20);
  nnbarStyle->SetMarkerSize(2);
  nnbarStyle->SetHistLineWidth(2.);
  nnbarStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  nnbarStyle->SetMarkerColor(kAzure+2);
  nnbarStyle->SetLineColor(kAzure+2);

  // get rid of X error bars (as recommended in ATLAS figure guidelines)
  nnbarStyle->SetErrorX(0.0001);
  // get rid of error bar caps
  nnbarStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  nnbarStyle->SetOptTitle(0);
  //nnbarStyle->SetOptStat(1111);
  nnbarStyle->SetOptStat(0);
  //nnbarStyle->SetOptFit(1111);
  nnbarStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  nnbarStyle->SetPadTickX(1);
  nnbarStyle->SetPadTickY(1);

  return nnbarStyle;




}  
