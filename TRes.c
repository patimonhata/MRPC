#include <fstream>
#include "ChMap.h" // homemade class

//function prototypes. Implementations after the main fuction.
std::array<std::vector<int>,2> setDiscriChOfRPC(const std::vector<string> allStrips, const std::vector<string>& discriInput, const ChMap& chMap);
std::array<double, 8> getTiming(std::vector<std::vector<double>>* ltdc, std::vector<std::vector<double>>* ttdc, int rch, int lch);
double calPercentile(double tot, std::vector<double> sortedTots);
TH2F* convert2IntegralvsTOF(TGraph *gTOTvsTOF, std::vector<int> decentEvents, std::vector<double> sortedTots);
TH2F* fillHistTotvsFraction(std::vector<double> sortedTots);
TH2F* subtractSlewing (TH2F* hIntegralvsTOF, TF1* f);

//main()
void TRes(int run, double fpeak_fix=0, double tdc_finetune=0)
{

  //E50 Carbonless-RPC
  //202401beamtest
  
  const double interval = 1.96625;
  const double fpeak = 409.313 + fpeak_fix;
  const double tdc_corr_fac = 1.04002 + tdc_finetune;



  std::vector<string> allStrips = {"1L","2L","3L","4L","5L","6L","7L","8L","9L","10L","11L","12L","13L","14L","15L","16L","1R","2R","3R","4R","5R","6R","7R","8R","9R","10R","11R","12R","13R","14R","15R","16R"};
  std::vector<int> watched_strips;
  int theStrip;
  ChMap chMap;
  const int n_of_discri_ch = 32;
  std::vector<string> discriInput;
  
  /*set a chMap*/
  if (1001 <= run && run <= 1018) {
    theStrip = 6; //the strip on which the trigger scintillators are.
    watched_strips = {5,6,7}; //neighbor strips that are monitored 
    discriInput = {"13L","","5R","","12L","","6R","","","6L","","12R","","5L","","13R","trig1","","trig2","","trig3","","trig4","","7L","14R","","7R","14L","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1019 <= run&& run <= 1040) {  
    theStrip = 13;
    watched_strips = {12,13,14};
    discriInput = {"13L","","5R","","12L","","6R","","","6L","","12R","","5L","","13R","trig1","","trig2","","trig3","","trig4","","7L","14R","","7R","14L","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1041 <= run && run <= 1075) {
    theStrip = 13;
    watched_strips = {12,13,14};
    discriInput = {"13L","","5R","14R","12L","","6R","14L","","6L","","12R","","5L","","13R","trig1","","trig2","","trig3","","trig4","","7L","","","7R","","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1076 <= run && run <= 1082){
    theStrip = 10;
    watched_strips = {9,10,11};
    discriInput = {"3L","","11R","2R","4L","","10R","2L","","10L","","4R","","11L","","3R","trig1","","trig2","","trig3","","trig4","","9L","","","9R","","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1083 <= run && run <= 1085){
    theStrip = 10;
    watched_strips = {9,10,11};
    discriInput = {"3L","","11R","2R","9L","","10R","2L","","10L","","9R","","11L","","3R","trig1","","trig2","","trig3","","trig4","","4L","","","4R","","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1086 <= run && run <= 1092) {  
    theStrip = 3;
    watched_strips = {2,3,4};
    discriInput = {"3L","","11R","2R","4L","","10R","2L","","10L","","4R","","11L","","3R","trig1","","trig2","","trig3","","trig4","","9L","","","9R","","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1093 <= run && run <= 1097){
    theStrip = 14;
    watched_strips = {13,14,15};
    discriInput = {"14L","","6R","15R","13L","","7R","15L","","7L","","13R","","6L","","14R","trig1","","trig2","","trig3","","trig4","","8L","","","8R","","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else if (1098 <= run && run <= 1099){
    theStrip = 7;
    watched_strips = {6,7,8};
    discriInput = {"14L","","6R","15R","8L","","7R","15L","","7L","","8R","","6L","","14R","trig1","","trig2","","trig3","","trig4","","13L","","","13R","","","RF","tag"};
    chMap.Set(discriInput, n_of_discri_ch);
  } else {
    // do nothing
  }
  
  int lch;
  int rch;
  int RFch;
  int tagch;
  //tigger ch
	const int nTrig = 4;
	int trigChs[nTrig];
	for( int i=0; i<nTrig; i++){ trigChs[i] = chMap.Get(Form("trig%d",i+1)); }
	//RFch
	RFch = chMap.Get("RF");
  tagch = chMap.Get("tag");
	//RPCch
	auto [chRPCL, chRPCR] = setDiscriChOfRPC(allStrips, discriInput, chMap);
  lch = chRPCL[theStrip - 1]; //vector is 0-start while the strip number is 1-start
  rch = chRPCR[theStrip - 1];

  const int Nbunch = 430;
  
  /*parameters for slewing correction*/
  const double loFitVal = 0;
  const double hiFitVal = 100;
  const double q = 100;

  /*tot cut parameter*/
  const double lowtot = 0;
  const double hightot = 40;
  
  double res[3], err;
  double res1[3], err1;
  
  bool decent_event_flag; //both sides of the strip is ringing & tot>0 & RF timing is whitin healthy range
  
  //inputfile path
  char filename[100];
  //sprintf(filename,"/home/daq1/work/RPC_202310/RPC1/Analyzer/root/run%04d.root",run);
  sprintf(filename,"root/run%04d.root",run);

  TFile *fin = new TFile(filename);
  TTree *tree = (TTree*)fin->Get("tree");

  std::vector<std::vector<double>> *ltdc = 0;
  std::vector<std::vector<double>> *ttdc = 0;
  TBranch *Bltdc = 0;
  TBranch *Bttdc = 0;
  tree->SetBranchAddress("ltdc", &ltdc, &Bltdc );
  tree->SetBranchAddress("ttdc", &ttdc, &Bttdc );

  std::vector<double> tots;
  std::vector<int> decentEvents; //the events in which both sides of the strip is ringing & RF timing is whitin healthy range

  TGraph *gtot_tof = new TGraph();

  TH2F *htot_fr = new TH2F(); 
  TH1F *htot = new TH1F("htot", "tot", 50, 0, 10);

  TH1F *htofRestored = new TH1F("htofRestored", "RFreal-RPCmean", 100, -2, 2);
  TH1F *htofBefore = new TH1F("htofBefore", "RF-RPCmean-fpeak", 10000, -200, 100);
  TH2F *htot_tof = new TH2F("htot_tof", "tot_RF-RPC", 800, -1, 10, 200, -1.5, 1.5);
  TH2F *htot_tof_f = new TH2F("htot_tof_f", "tot_RF-RPC-f", 200, -1, 10, 300, -1.5, 1.5);     

  //TH2F *hfr_tof = new TH2F("hfr_tof", "fraction_RF-RPC", 53, -3, 103, 300, -1.5, 1.5);
  TH2F *h3 = new TH2F();

  TH1D *htof1 = new TH1D("htof1", "RFreal-RPC-slewing", 100, -2, 2);  
  TH2F *hfr_tof1 = new TH2F(); 

  //TH1F *h01 = new TH1F("h01", "trig2-trig1", 200, -100, 100);
  //TH1F *h02 = new TH1F("h02", "trig3-trig1", 200, -100, 100);
  //TH1F *h03 = new TH1F("h03", "trig3-trig2", 200, -100, 100);
  
  double tot, totL, totR;
  double leadL, leadR;
  double trailL, trailR;
  double leadmean;
  double rfRecorded;
  double rfReal;
  double trig1;
  double trig2;
  double trig3;
  double trig4;
  double eff,  effe;
  double effR, effLe;
  double effL, effRe;
  int ev = 0; //triggers and RF captured
  int Rev = 0; //at least right side ringing
  int Lev = 0; //at least left side ringing
  int rpcev = 0; //both side ringing
  
  double Fraction;
  int k = 0;
  
  const int n = tree->GetEntries(); //total events

  ///////////////1////////////////
  /* Efficiency evaluation */
  for(int i=0; i<n; i++){

    tree->GetEntry(i);

    //ensure that trigger scinti's & RF are captured
    if( ltdc->at(trigChs[0]).size() < 1 || ttdc->at(trigChs[0]).size() < 1 ) continue;
    if( ltdc->at(trigChs[1]).size() < 1 || ttdc->at(trigChs[1]).size() < 1 ) continue;
    if( ltdc->at(trigChs[2]).size() < 1 || ttdc->at(trigChs[2]).size() < 1 ) continue;
    if( ltdc->at(trigChs[3]).size() < 1 || ttdc->at(trigChs[3]).size() < 1 ) continue;
    if(ltdc->at(RFch).size()<1 || ttdc->at(RFch).size()<1)continue;
    
    ev +=1;  //trigger event

    //rpc
    if(ltdc->at(lch).size()>0 && ttdc->at(lch).size()>0) Lev+=1;
    if(ltdc->at(rch).size()>0 && ttdc->at(rch).size()>0) Rev+=1;
    if(ltdc->at(lch).size()>0 && ttdc->at(lch).size()>0
       && ltdc->at(rch).size()>0 && ttdc->at(rch).size()>0) rpcev +=1;
    
    // maybe we shuould increment rpcev when the event is decent (~10 counts difference) 

  }


  ///////////////2////////////////
  // fill htot, htot_tof, htofRestored.
  for(int i=0; i<n; i++){

    tree->GetEntry(i);
    
    decent_event_flag = false;
    
    //ensure that trigger scinti's & RF are captured
    if( ltdc->at(trigChs[0]).size() < 1 || ttdc->at(trigChs[0]).size() < 1 ) continue;
    if( ltdc->at(trigChs[1]).size() < 1 || ttdc->at(trigChs[1]).size() < 1 ) continue;
    if( ltdc->at(trigChs[2]).size() < 1 || ttdc->at(trigChs[2]).size() < 1 ) continue;
    if( ltdc->at(trigChs[3]).size() < 1 || ttdc->at(trigChs[3]).size() < 1 ) continue;
    if(ltdc->at(RFch).size()<1 || ttdc->at(RFch).size()<1)continue;
    //tagger cut
    //    if( ltdc->at(tagch).size() < 1 || ttdc->at(tagch).size() < 1 ) continue;
    //    if( ltdc->at(tagch).at(0) < 205 || ttdc->at(tagch).at(0)() >225 ) continue;
 
    //ensure that both leading and trailing of both side of the strip are captured
    if(ltdc->at(lch).size()<1 || ttdc->at(lch).size()<1) continue;
    if(ltdc->at(rch).size()<1 || ttdc->at(rch).size()<1) continue;

    // insert the signal timings or tots into variables
    auto [leadL, trailL, totL, leadR, trailR, totR, leadmean, tot] = getTiming(ltdc, ttdc, rch, lch);
    rfRecorded = ltdc->at(RFch).at(0) * tdc_corr_fac;
    
    if(tot<0) continue;

    htofBefore->Fill(rfRecorded-leadmean-fpeak);

    // convert a rfRecorded to the real RF timing  
    for(int j=0; j<Nbunch; j++){
      if( fabs(rfRecorded - leadmean - fpeak - interval*j) <= interval/2 ){
        rfReal = rfRecorded - fpeak - interval*j;
	      decent_event_flag = true;
      }
    }

    if(decent_event_flag == false) continue;
    decentEvents.push_back(i);
    tots.push_back(tot);

    htot->Fill(tot);
    htot_tof->Fill(tot,rfReal-leadmean);
    htofRestored->Fill(rfReal-leadmean);
    gtot_tof->SetPoint(i,tot,rfReal-leadmean);
  }//end of eventloop

  //cout << decentEvents.size() << endl;
  sort(tots.begin(),tots.end());

  
  TCanvas *c0 = new TCanvas("c0","c0",3,2,1400,700);
  c0->Divide(5,2);

  
  c0->cd(1);
  htot->Draw("box");

  c0->cd(2);
  htofRestored->Draw("box");
  TF1 *Gaus0 = new TF1("Gaus0","gaus(0)");
  htofRestored->Fit("Gaus0","","", -2, 2);
  Gaus0->GetParameters(res);
  err = Gaus0->GetParError(2);
 
  c0->cd(3);
  htot_tof->Draw("colz");

  c0->cd(10);
  htofBefore->Draw("box");
  
  h3 = convert2IntegralvsTOF(gtot_tof, decentEvents, tots);/////////////my convert///////////////////
  c0->cd(4);
  h3->Draw("colz");

  TH1D* hh3_pfx = (TH1D*)h3->ProfileX("hh3_pfx");
  hh3_pfx->SetTitle("slew. in low<TOT<q");
  TH1D* hh4_pfx = (TH1D*)h3->ProfileX("hh4_pfx");
  hh4_pfx->SetTitle("slew. in q<TOT<high");
  c0->cd(5);
  hh3_pfx->Draw();
  TF1 *pol=new TF1("pol","[0]+[1]*x+[2]*x*x",0,700);
  hh3_pfx->Fit("pol","","", loFitVal, q);
  c0->cd(6);
  hh4_pfx->Draw();
  TF1 *pol2=new TF1("pol2","[0]+[1]*x+[2]*x*x",0,700);
  hh4_pfx->Fit("pol2","","", q, hiFitVal);



  hfr_tof1 = subtractSlewing(h3, pol);/////////////////my subtractSlewing///////////////////////

  htof1 = dynamic_cast<TH1D*>(hfr_tof1->ProjectionY("htof1")->Rebin(4));

  c0->cd(7);
  hfr_tof1->Draw("colz");
  //htot_tof_f->Draw("colz");


  c0->cd(8);
  htof1->Draw("box");
  TF1 *Gaus1 = new TF1("Gaus1","gaus(0)");
  htof1->Fit("Gaus1","","", -2, 2);
  Gaus1->GetParameters(res1);
  err1 = Gaus1->GetParError(2);

  c0->cd(9);
  htot_fr = fillHistTotvsFraction(tots); /////////////////my fillHistTotvsFraction///////////////////////
  htot_fr->Draw("box");


  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.16);


  /// last ///

    
   
    cout << "--------------------------------------------------" << endl;
    cout << "total event = " << n << endl;
    cout << "trig event = " << ev << endl;
    cout << "Left event = " << Lev << endl;
    cout << "Right event = " << Rev << endl;
    cout << "rpc event = " << rpcev << endl;

    effL = (double)Lev/ev *100;
    effLe = sqrt(effL/100*(1-effL/100))/sqrt(ev)*100;
    cout << "L effciency = " << effL << " pm " << effLe << " %" << endl;
    effR = (double)Rev/ev *100;
    effRe = sqrt(effR/100*(1-effR/100))/sqrt(ev)*100;
    cout << "R effciency = " << effR << " pm " << effRe << " %" << endl;
    eff = (double)rpcev/ev *100;
    effe = sqrt(eff/100*(1-eff/100))/sqrt(ev)*100;
    cout << "effciency = " << eff << " pm " << effe << " %" << endl;

    cout << "Resolution (before correction) = " << res[2]*1000 << " pm " << err*1000 << " ps" << endl;
    cout << "Resolution (after correction) = " << res1[2]*1000 << " pm " << err1*1000 << " ps" << endl;

    //htot_tof->SetStats(0);    
    //h3->SetStats(0);
    //hfr_tof1->SetStats(0);
    //htof1->SetStats(0);

    /*
    char outfile[100];
    sprintf(outfile,"result/run%04d.txt",run);
    ofstream outputfile(outfile);
    outputfile << "-------------------------------------------------- \n"
    << "total event = " << n << "\n"
    << "trig event = " << ev << "\n"
    << "Left event = " << Lev << "\n"
    << "Right event = " << Rev << "\n"
    << "rpc event = " << rpcev << "\n";

    effL = (double)Lev/ev *100;
    effLe = sqrt(effL/100*(1-effL/100))/sqrt(ev)*100;
    outputfile << "L effciency = " << effL << " pm " << effLe << " %" << "\n";
    effR = (double)Rev/ev *100;
    effRe = sqrt(effR/100*(1-effR/100))/sqrt(ev)*100;
    outputfile << "R effciency = " << effR << " pm " << effRe << " %" << "\n";
    eff = (double)rpcev/ev *100;
    effe = sqrt(eff/100*(1-eff/100))/sqrt(ev)*100;
    outputfile << "effciency = " << eff << " pm " << effe << " %" << "\n";

    outputfile << "Resolution (before correction) = " << res[2]*1000 << " pm " << err*1000 << " ps" << "\n";
    outputfile << "Resolution (after correction) = " << res1[2]*1000 << " pm " << err1*1000 << " ps" << "\n";;
    outputfile << "-------------------------------------------------- \n";

    outputfile.close();
    
    char outname[100];
    sprintf(outname, "pic/run%04d.jpg",run);
    c0->Print(outname);
    */

}// end of main()











std::array<std::vector<int>,2> setDiscriChOfRPC(const std::vector<string> allStrips, const std::vector<string>& discriInput, const ChMap& chMap) {
    std::vector<int> chRPCL;
    std::vector<int> chRPCR;

    for (string s : allStrips){
        bool inputted = std::find(discriInput.begin(), discriInput.end(), s) != discriInput.end(); //discriInput の中に s が存在するか
        if (inputted){
        if (s.back() == 'L'){
            chRPCL.push_back(chMap.Get(s)); //ストリップsが discri に入力されていれば ch を返す
        } else {
            chRPCR.push_back(chMap.Get(s));
        }
        } else {
        if (s.back() == 'L'){
            chRPCL.push_back(chMap.Get("")); //ストリップsが discri に入力されていれば ch を返す
          } else {
            chRPCR.push_back(chMap.Get(""));
        }
        // これでは discri 全chに挿していて、かつ使っていないストリップがあるときにエラーで終了してしまう。改善の必要あり。
        }
    }

    return {std::move(chRPCL), std::move(chRPCR)}; //{chRPCL, chRPCR} で初期化
};

std::array<double, 8> getTiming(std::vector<std::vector<double>>* ltdc, std::vector<std::vector<double>>* ttdc, int rch, int lch){
    // use like [leadL, trailL, totL, leadR, trailR, totR, leadmean, tot] = getTiming(tree);

    double tdc_corr_fac = 1.04002;
    
    double leadL = ltdc->at(lch).at(0) * tdc_corr_fac;
    double trailL = ttdc->at(lch).at(0) * tdc_corr_fac;
    double totL = leadL - trailL;
    double leadR = ltdc->at(rch).at(0) * tdc_corr_fac;
    double trailR = ttdc->at(rch).at(0) * tdc_corr_fac;
    double totR = leadR-trailR; 
    double leadmean = (leadL + leadR) /2;
    double tot = (totL + totR) /2;

    return {
        std::move(leadL),
        std::move(trailL), 
        std::move(totL), 
        std::move(leadR), 
        std::move(trailR), 
        std::move(totR),
        std::move(leadmean),
        std::move(tot)
    };
}

double calPercentile(double tot, std::vector<double> sortedTots){
    auto it = std::lower_bound(sortedTots.begin(), sortedTots.end(), tot); // get iterator
    double percentile = (it - sortedTots.begin()) * 100.0 / sortedTots.size(); // calculate percentile of tot in the tots
    return percentile;
}

TH2F* convert2IntegralvsTOF(TGraph *gTOTvsTOF, std::vector<int> decentEvents, std::vector<double> sortedTots){
  TH2F* hIntegralvsTOF = new TH2F("h3","percentile_TOF", 106,-3,103, 200,-1.5,1.5);
  double tot;
  double percentile;

  //sort(tots.begin(), tots.end());

  for (int i : decentEvents) {
    tot = gTOTvsTOF->GetPointX(i);
    percentile = calPercentile(tot, sortedTots);
    hIntegralvsTOF->Fill(percentile, gTOTvsTOF->GetPointY(i));
  }
  
  return std::move(hIntegralvsTOF);
}

TH2F* fillHistTotvsFraction(std::vector<double> sortedTots){
  TH2F* hTOTvsFraction = new TH2F("htot_fr", "tot_percentile", 400, -1, 40, 10000, 0, 100);
  double tot;
  double percentile;

  //sort(tots.begin(), tots.end());

  for (double tot : sortedTots) {
    percentile = calPercentile(tot, sortedTots);
    hTOTvsFraction->Fill(tot, percentile);
  }

  return std::move(hTOTvsFraction);
}

//for now, does not support slewing correction which subtract different function in different sections of tot. 
//(dividing tot into 2 section and subtract different function in each section is supposed)
TH2F* subtractSlewing (TH2F* hIntegralvsTOF, TF1* f){  
    int nBinsX = hIntegralvsTOF->GetNbinsX();
    int nBinsY = hIntegralvsTOF->GetNbinsY();
    TH2F* hIntegralvsTOF_subtracted = new TH2F("hIntegralvsTOF_subtracted","after slew. corr.", 100,0,100,
                                    nBinsY,
                                    hIntegralvsTOF->GetYaxis()->GetBinLowEdge(1),
                                    hIntegralvsTOF->GetYaxis()->GetBinLowEdge(nBinsY + 1)
                                    ); //create a TH2F that have same bins

    for (int i=0; i<nBinsX; i++){ //i,j = the index of x, y, respectively;
        for (int j=0; j<nBinsY; j++){
            int bin = hIntegralvsTOF->GetBin(i,j); //get global bin number from xbin and ybin
            int binContent = hIntegralvsTOF->GetBinContent(bin); //the number of entries in that bin
            if (binContent ==0) continue;

            double x_Integral = hIntegralvsTOF->GetXaxis()->GetBinCenter(i); //the center of the i-th bin in x direction;
            double y = hIntegralvsTOF->GetYaxis()->GetBinCenter(j); //the center of the j-th bin in y direction;

            //fill for the binContent times
            for (int a=0; a<binContent; a++){
                hIntegralvsTOF_subtracted->Fill(x_Integral, y - f->Eval(x_Integral));
            }

        }
    }

    return std::move(hIntegralvsTOF_subtracted);
}