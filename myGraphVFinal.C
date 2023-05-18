void format_line(TAttLine* line,int col,int sty){ // Función de estilo
    line->SetLineWidth(1.0); line->SetLineColor(col);
    line->SetLineStyle(sty);}

double the_gauspparerf(double* vars, double* pars){
    return pars[0]*TMath::Gaus(vars[0],pars[1],pars[2])+
        pars[3]+pars[4]*vars[0]+pars[5]*vars[0]*vars[0]+
        pars[6]*TMath::Erf((vars[0]-pars[7])/pars[8]);}

void myGraphVFinal(){

    gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111); gStyle->SetStatBorderSize(0);
    // Open the input file containing the TTree
    TFile* infile = TFile::Open("reducetree_Bujk_AOD_3_best1.root");

    // Get the TTree from the input file
    TTree* tree = (TTree*)infile->Get("butree");
    
    // Set the branch
    double B_mass;
    tree->SetBranchAddress("B_mass", &B_mass);

    // Create a Histogram
    int num_bins = 100;
    float min_value = 5;
    float max_value = 5.6;
    auto hist = new TH1F("Histogram", "Histogram Title", num_bins, min_value, max_value);

    // Fill the Hitogram with the values from the TTree branches
    // tree->GetEntries(): 770929
    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i); // Output: 132
        hist->Fill(B_mass);
        //cout << "Entries: " << h->GetEntries() << endl;
        //cout << "B_mass value: " << B_mass << endl;
    }

    // Create a canvas and draw the TGraph

    TCanvas* canvas = new TCanvas("canvas", "Graph", 800, 600);

    //hist->Draw();
    // Customize histogram appearance
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->GetXaxis()->SetNdivisions(5, 5, 0);
    hist->GetXaxis()->SetRangeUser(5.0, 5.6);

    hist->SetMinimum(1e2); // The minimum value should be greater than zero
    hist->SetMaximum(1e5); 
    

    canvas->SetLogy(true);

    auto parabola = new TF1("parabola","[0]+[1]*x+[2]*x**2",0,20);
    format_line(parabola,kBlack,0);

    auto gaussian = new TF1("gaussian","[0]*TMath::Gaus(x,[1],[2])",0,20);
    format_line(gaussian,kRed,0);

    auto erf = new TF1("erf", "[0]*TMath::Erf((x-[1])/[2])", 0, 20);
    format_line(erf,kYellow,1);

    auto gauspparerf = new TF1("gausppar",the_gauspparerf,5.0,5.6,9);// El 9 es el número de parámetros
    format_line(gauspparerf,kBlue,0);

    gauspparerf->SetParNames("Norm","Mean","Sigma","a","b","c","Norm_erf","Mean_erf","Sigma_erf"); // Nombramos los parámetros

    gauspparerf->SetParameter(0,5.47e4);
    gauspparerf->SetParameter(1,5.3);
    gauspparerf->SetParameter(2,1.8e-2);
    int npar=gauspparerf->GetNpar();
    for (int ipar=3;ipar<npar-3;++ipar) gauspparerf->SetParameter(ipar,1);
    gauspparerf->SetParameter(6,1e4);
    gauspparerf->SetParameter(7,5.3);
    gauspparerf->SetParameter(8,0.1);
    
    hist->Fit(gauspparerf, "R");
    gauspparerf->DrawClone("Same");

    // Given the parameters to the gaussian and parabola functions
    for (int ipar=0;ipar<3;ipar++){
        gaussian->SetParameter(ipar,gauspparerf->GetParameter(ipar));
        parabola->SetParameter(ipar,gauspparerf->GetParameter(ipar+3));
        erf->SetParameter(ipar,gauspparerf->GetParameter(ipar+6));
    }

    hist->DrawClone("PE");
    parabola->DrawClone("Same");
    gaussian->DrawClone("Same");
    erf->DrawClone("Same");
    
    canvas->SaveAs("histogram.png");

    // Close the input file
    //infile->Close();
}