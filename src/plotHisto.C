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


  // set label offset
  nnbarStyle->SetLabelOffset(0.01,"xyz");

  // use large fonts
  Int_t font=62; //  helvetica-medium-r-normal "Arial"
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




  //gROOT->Reset();
  gROOT->SetStyle("NNbar");
  gROOT->ForceStyle("NNbar");
 
//......................................................

  // Draw histos filled by Geant4 simulation 
  //  

  // Open file filled by Geant4 simulation 
  TFile f("output.root");

  // Create a canvas and divide it into 5x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 200, 200, 1000, 1000);
  c1->UseCurrentStyle();

  c1->Divide(2,5);
 

 
  // Draw histogram in pad 1
  c1->cd(1);
  TH1D* hist1 = (TH1D*)f.Get("scintphotons_A");
  hist1->Draw("HIST");
  
  // Draw histogram in pad 2
  c1->cd(2);
  TH1D* hist2 = (TH1D*)f.Get("scintphotons_B");
  hist2->Draw("HIST");
  
  // Draw histogram in pad 3
  TH1D* hist3 = (TH1D*)f.Get("scintphotons_C");
  c1->cd(3);
  hist3->Draw("HIST");
  
  // Draw histogram in pad 4
  c1->cd(4);
  TH1D* hist4 = (TH1D*)f.Get("scintphotons_D");
  hist4->Draw("HIST");

  // Draw histogram in pad 5
  c1->cd(5);
  TH1D* hist5 = (TH1D*)f.Get("scintphotons_E");
  hist5->Draw("HIST");
  
  // Draw histogram in pad 6
  c1->cd(6);
  TH1D* hist6 = (TH1D*)f.Get("scintphotons_F");
  hist6->Draw("HIST");
  
  // Draw histogram in pad 7
  TH1D* hist7 = (TH1D*)f.Get("scintphotons_G");
  c1->cd(7);
  hist7->Draw("HIST");
  
  // Draw histogram in pad 8
  c1->cd(8);
  TH1D* hist8 = (TH1D*)f.Get("scintphotons_H");
  hist8->Draw("HIST");

  // Draw histogram in pad 9
  TH1D* hist9 = (TH1D*)f.Get("scintphotons_I");
  c1->cd(9);
  hist9->Draw("HIST");
  
  // Draw histogram in pad 10
  c1->cd(10);
  TH1D* hist10 = (TH1D*)f.Get("scintphotons_J");
  hist10->Draw("HIST");






  
  c1->SaveAs("name.png");

  return nnbarStyle;




}  
