// ROOT macro file for plotting example B4 histograms 
// 
// Can be run from ROOT session:
// root[0] .x plotHisto.C

// Make filename user input

{

//................
// Define Style
//
TStyle *nnbarStyle = new TStyle("NNbar","NNbar style");

// use plain black on white colors
Int_t icol=0; 
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
Int_t font=62; 
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
//nnbarStyle->SetMarkerStyle(20);
//nnbarStyle->SetMarkerSize(2);
//nnbarStyle->SetHistLineWidth(2.);
//nnbarStyle->SetLineStyleString(2,"[12 12]"); 
//nnbarStyle->SetMarkerColor(kAzure+2);
//nnbarStyle->SetLineColor(kAzure+2);

// get rid of X error bars 
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

gROOT->SetStyle("NNbar");
gROOT->ForceStyle();


//......................................................

// Draw histos filled by Geant4 simulation 

///// Test
filename = "../output/calo-sim.root";
outfile_a = "pres-1D.png";

// Open file filled by Geant4 simulation 
TFile f(filename);

// Create a canvas and divide it into 2x2 pads
//TCanvas *c1 = new TCanvas("c1", "", 585, 276, 961, 646);
TCanvas *c1 = new TCanvas("c1", "", 500, 500, 1300, 900);
c1->Divide(1,2);

////// Draw c1 plots
// 1    2
// 3    4
// 5    6

TH1D* h1 = (TH1D*)f.Get("EdepScint");
TH1D* h2 = (TH1D*)f.Get("EdepAbs");


c1->cd(1);//->SetRightMargin(0.18);
//gPad->SetLogy(1);
h1->SetTitle("");
h1->GetXaxis()->SetTitle("Energy Deposited in Scint [MeV]");
h1->GetYaxis()->SetTitle("Num Events");
//h1->SetMaximum(2e4);
//h1->SetMinimum(1e-1);
h1->SetLineColor(kBlack);
h1->SetMarkerSize(.5);
h1->SetMarkerStyle(20);
h1->Draw("EHIST");


c1->cd(2);//->SetRightMargin(0.18);
//gPad->SetLogy(1);
h2->SetTitle("");
h2->GetXaxis()->SetTitle("Energy Deposited in Absorber [MeV]");
h2->GetYaxis()->SetTitle("Num Events");
//h2->SetMaximum(2e4);
//h2->SetMinimum(1e-1);
h2->SetLineColor(kBlack);
h2->SetMarkerSize(.5);
h2->SetMarkerStyle(20);
h2->Draw("EHIST");

c1->Print(outfile_a);
}  
