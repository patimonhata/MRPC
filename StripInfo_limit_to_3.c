#include <algorithm>
#include "TStyle.h"
#include "TCanvas.h"
#include "stdbool.h"

#include "ChMap.h"

//function prototypes
void setStyle();
std::array<std::vector<int>,2> setDiscriChOfRPC(const std::vector<string> allStrips, const std::vector<string>& discriInput, const ChMap& chMap);

//main function
void StripInfo_limit_to_3(int run){
	setStyle();

  	std::vector<string> allStrips = {"1L","2L","3L","4L","5L","6L","7L","8L","9L","10L","11L","12L","13L","14L","15L","16L","1R","2R","3R","4R","5R","6R","7R","8R","9R","10R","11R","12R","13R","14R","15R","16R"};
  	std::vector<int> watched_strips;
	ChMap chMap;
  	const int n_of_discri_ch = 32;
  	std::vector<string> discriInput;

  	if (1001 <= run && run <= 1018) {
    	discriInput = {"13L","","5R","","12L","","6R","","","6L","","12R","","5L","","13R","trig1","","trig2","","trig3","","trig4","","7L","14R","","7R","14L","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {5, 6, 7};
  	} else if (1019 <= run&& run <= 1040) {  
    	discriInput = {"13L","","5R","","12L","","6R","","","6L","","12R","","5L","","13R","trig1","","trig2","","trig3","","trig4","","7L","14R","","7R","14L","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {12, 13, 14};
  	} else if (1041 <= run && run <= 1075) {
    	discriInput = {"13L","","5R","14R","12L","","6R","14L","","6L","","12R","","5L","","13R","trig1","","trig2","","trig3","","trig4","","7L","","","7R","","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {12,13,14};
  	} else if (1076 <= run && run <= 1082){
    	discriInput = {"3L","","11R","2R","4L","","10R","2L","","10L","","4R","","11L","","3R","trig1","","trig2","","trig3","","trig4","","9L","","","9R","","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {9,10,11};
  	} else if (1083 <= run && run <= 1085){
    	discriInput = {"3L","","11R","2R","9L","","10R","2L","","10L","","9R","","11L","","3R","trig1","","trig2","","trig3","","trig4","","4L","","","4R","","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {9,10,11};
  	} else if (1086 <= run && run <= 1092) {  
    	discriInput = {"3L","","11R","2R","4L","","10R","2L","","10L","","4R","","11L","","3R","trig1","","trig2","","trig3","","trig4","","9L","","","9R","","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {2,3,4};
  	} else if (1093 <= run && run <= 1097){
    	discriInput = {"14L","","6R","15R","13L","","7R","15L","","7L","","13R","","6L","","14R","trig1","","trig2","","trig3","","trig4","","8L","","","8R","","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {13,14,15};
  	} else if (1098 <= run && run <= 1099){
    	discriInput = {"14L","","6R","15R","8L","","7R","15L","","7L","","8R","","6L","","14R","trig1","","trig2","","trig3","","trig4","","13L","","","13R","","","RF","tag"};
    	chMap.Set(discriInput, n_of_discri_ch);
    	watched_strips = {6,7,8};
  	} else {
    	// do nothing
  	}


	//tigger ch
	int nTrig = 4;
	int trigChs[nTrig];
	for( int i=0; i<nTrig; i++){ trigChs[i] = chMap.Get(Form("trig%d",i+1)); }
	//RFch
	int RFch = chMap.Get("RF");
	//RPCch
	auto [chRPCL, chRPCR] = setDiscriChOfRPC(allStrips, discriInput, chMap);
	//for (int i : chRPCL) cout << "chRPCL : "<< i << endl;
	//for (int i : chRPCR) cout << "chRPCR : "<< i << endl;

	//tdc cut
	double tdcrangelowL = 200;
	double tdcrangehighL = 300;
	double tdcrangelowR = 200;
	double tdcrangehighR = 300;

	TFile *fin = new TFile(Form("root/run%04d.root",run));
	cout << Form("root/run%04d",run) << " was opened" << endl;
	TTree *tree = (TTree*)fin->Get("tree");

	std::vector<std::vector<double>> *ltdc = 0;
	std::vector<std::vector<double>> *ttdc = 0;
	TBranch *Bltdc = 0;
	TBranch *Bttdc = 0;
	tree->SetBranchAddress("ltdc", &ltdc, &Bltdc );
	tree->SetBranchAddress("ttdc", &ttdc, &Bttdc );


	//definitions of histograms
	TH1F *hHitChL = new TH1F("hHitChL","L Hit Map; Ch No.",16,0,16);
	TH1F *hHitChR = new TH1F("hHitChR","R Hit Map; Ch No.",16,0,16);
	TH2F *hCh_LvsR = new TH2F("Correlation of Left strip and Right strip","; L Ch No.; R Ch No.",16,0,16,16,0,16);//16,0,16,16,0,16);
	TH2F *hCh_L1vsL2 = new TH2F("Correlation of Left strips", "Correlation of Left strips; first strip ; second strip",16,0,16,16,0,16);
	TH2F *hCh_R1vsR2 = new TH2F("Correlation of Right strips", "Correlation of Right strips; first strip ; second strip",16,0,16,16,0,16);

	////leading histograms
	TH1F *hltdcL[chRPCL.size()];
	for( int i=0; i<chRPCL.size(); i++){
		hltdcL[i] = new TH1F(Form("hltdcL%d",i),Form("ltdc L%d; ltdc [ns]",i+1),100,200,300);
	}
  
	TH1F *hltdcR[chRPCR.size()];
	for( int i=0; i<chRPCR.size(); i++){
    	hltdcR[i] = new TH1F(Form("hltdcR%d",i),Form("ltdc R%d; ltdc [ns]",i+1),100,200,300);
  	}

  ////TOT histograms
  TH1F *htotL[chRPCL.size()];
  for( int i=0; i<chRPCL.size(); i++){
    htotL[i] = new TH1F(Form("htotL%d",i),Form("TOT L%d; TOT [ns]",i+1),50,-1,10);
  }

  TH1F *htotR[chRPCR.size()];
  for( int i=0; i<chRPCR.size(); i++){
    htotR[i] = new TH1F(Form("htotR%d",i),Form("TOT R%d; TOT [ns]",i+1),50,-1,10);
  }

  ////ltdc of trigger
  TH1F *hltdc_trig[nTrig]; 
  for( int i=0; i<nTrig; i++){
    hltdc_trig[i]= new TH1F(Form("hltdc_trig[%d]",i), Form("ltdc_trig%d; ltdc [ns]",i+1),50,215,230);
  }

  //// hist of multiplicity
  TH1F *hMltptyL = new TH1F("hist of multiplicity L","hist of multiplicity L",4,-0.5,3.5);
  TH1F *hMltptyR = new TH1F("hist of multiplicity R","hist of multiplicity R",4,-0.5,3.5);

  //definitions of variables
  double totL;
  double totR;
  int ChMaxtotL, ChMaxtotR;
  double MaxtotL, MaxtotR;
  double leadL, leadR, leadTrig;
  double trailL, trailR, trailTrig;
  bool eventflagL, eventflagR;
  bool evflagL, evflagR;
  int ev=0, evL=0, evR=0;
  int sevL = 0, sevR = 0;
  int trev=0;

  
  const int n = tree->GetEntries();
  //const int n = 10000;

	//event loop
  for(int i=0; i<n; i++){
    if(i%10000 == 0) cout << "event num = " << i << endl;
    tree->GetEntry(i);

    // cout << i << endl;
    for (int i=0; i<nTrig; i++){
      if( ltdc->at(trigChs[i]).size()<1 || ttdc->at(trigChs[i]).size()<1 ) continue;
    }
    if ( ltdc->at(RFch).size()<1 || ttdc->at(RFch).size()<1 ) continue; //efficinecy に影響なかった

    //if( ltdc->at(TagCh).size()<1 || ttdc->at(TagCh).size()<1 ) continue;
    //if( ltdc->at(TagCh).at(0)<80 || ltdc->at(TagCh).at(0)>140 ) continue;


  //cout << i << endl;
    std::vector<int> HitChL;
    std::vector<int> HitChR;

    MaxtotL = 0;
    MaxtotR = 0;
    ChMaxtotL = 100;
    ChMaxtotR = 100;
    eventflagL = false;
    eventflagR = false;
    evflagL=false;
    evflagR=false;
    trev++;
      
    for(int j=0; j<chRPCL.size(); j++){
      if( ltdc->at(chRPCL[j]).size() < 1 || ttdc->at(chRPCL[j]).size() < 1 ) continue;
      leadL = 100;
      trailL = 100;
      leadL = ltdc->at(chRPCL[j]).at(0);
      trailL = ttdc->at(chRPCL[j]).at(0);
      totL = leadL - trailL;
      hltdcL[j]->Fill(leadL);
      if( totL<0 ) continue;
      //cout << "totL" << trailL << endl;
      htotL[j]->Fill(totL);
      //HitChL.push_back(j);
      
      for (int ss : watched_strips){
        if(j+1==ss && leadL>tdcrangelowL && leadL<tdcrangehighL){ 
          eventflagL=true;
          HitChL.push_back(j);
        }
      }
      if(j+1==watched_strips[1]&&leadL>tdcrangelowL&&leadL<tdcrangehighL) evflagL=true;
      if( totL > MaxtotL){
	      MaxtotL = totL;
	      ChMaxtotL = j;
      }
    }

    for(int j=0; j<chRPCR.size(); j++){
      if( ltdc->at(chRPCR[j]).size() < 1 || ttdc->at(chRPCR[j]).size() < 1 ) continue;
      
      leadR = 100;
      trailR = 100;
      leadR = ltdc->at(chRPCR[j]).at(0);
      trailR = ttdc->at(chRPCR[j]).at(0);
      totR = leadR-trailR;
      hltdcR[j]->Fill(leadR);

      htotR[j]->Fill(totR);
      //HitChR.push_back(j);
      
      for (int ss : watched_strips){
        if(j+1==ss && leadR>tdcrangelowR && leadR<tdcrangehighR){ 
          eventflagR=true;
          HitChR.push_back(j);
        }
      }
      if(j+1==watched_strips[1] && leadR>tdcrangelowR && leadR<tdcrangehighR) evflagR=true;
      if( totR > MaxtotR){
	      MaxtotR = totR;
	      ChMaxtotR = j;
      }
    }

    for(int j=0; j<HitChL.size(); j++){
      hHitChL->Fill(HitChL[j]);
    }
    for(int j=0; j<HitChR.size(); j++){
      hHitChR->Fill(HitChR[j]);
    }

    ////fill hltdc_trig
    for (int i=0; i<nTrig; i++){
      if( ltdc->at(trigChs[i]).size()  > 0 && ttdc->at(trigChs[i]).size() > 0 ){
        leadTrig = 100;
        trailTrig = 100;
        leadTrig = ltdc->at(trigChs[i]).at(0);
        trailTrig = ttdc->at(trigChs[i]).at(0);
        hltdc_trig[i]->Fill(leadTrig);//-ltdc->at(trigChs[1]).at(0)); /////////////////////////////
      };
    }

    //// fill hist if multiplicity
    hMltptyL->Fill(HitChL.size()-0.5);
    hMltptyR->Fill(HitChR.size()-0.5);

    if(eventflagL) evL++;
    if(eventflagR) evR++;
    if(eventflagL && eventflagR) ev++;
    if(evflagL) sevL++;
    if(evflagR) sevR++;

    hCh_LvsR->Fill(ChMaxtotL, ChMaxtotR);
    if(HitChL.size()==2){
      hCh_L1vsL2->Fill(HitChL[0], HitChL[1]);
    }
    if(HitChR.size()==2){
      hCh_R1vsR2->Fill(HitChR[0], HitChR[1]);
    }
  }





  // ltdc left canvas
  TCanvas *c0 = new TCanvas("c0","ltdc distributions on Left",3,2,1200,900);
  c0->Divide(4,4);  
  for( int i=0; i<chRPCL.size(); i++){
    c0->cd(i+1);
    hltdcL[i]->Draw("");
    hltdcL[i]->GetYaxis()->SetLabelSize(0.05);
    hltdcL[i]->GetXaxis()->SetLabelSize(0.07);
    hltdcL[i]->GetXaxis()->SetTitleSize(0.06);
    hltdcL[i]->GetXaxis()->SetTitleOffset(0.9);
    hltdcL[i]->GetXaxis()->SetNdivisions(505);
    hltdcL[i]->GetYaxis()->SetNdivisions(505);
  }

  // ltdc right canvas
  TCanvas *c1 = new TCanvas("c1","ltdc distributions on Right",3,2,1200,900);
  c1->Divide(4,4);  
  for( int i=0; i<chRPCR.size(); i++){
    c1->cd(i+1);
    hltdcR[i]->Draw("");
    hltdcR[i]->GetYaxis()->SetLabelSize(0.05);
    hltdcR[i]->GetXaxis()->SetLabelSize(0.07);
    hltdcR[i]->GetXaxis()->SetTitleSize(0.06);
    hltdcR[i]->GetXaxis()->SetTitleOffset(0.9);
    hltdcR[i]->GetXaxis()->SetNdivisions(505);
    hltdcR[i]->GetYaxis()->SetNdivisions(505);
  }

  // tot left canvas
  TCanvas *c2 = new TCanvas("c2","tot distribution on Left",3,2,1200,900);
  c2->Divide(4,4);  
  for( int i=0; i<chRPCL.size(); i++){
    c2->cd(i+1);
    htotL[i]->Draw("");
    htotL[i]->GetYaxis()->SetLabelSize(0.05);
    htotL[i]->GetXaxis()->SetLabelSize(0.07);
    htotL[i]->GetXaxis()->SetTitleSize(0.06);
    htotL[i]->GetXaxis()->SetTitleOffset(0.9);
    hltdcL[i]->GetXaxis()->SetNdivisions(505);
    hltdcL[i]->GetYaxis()->SetNdivisions(505);
  }


  // tot right canvas
  TCanvas *c3 = new TCanvas("c3","tot distributioin on Right",3,2,1200,900);
  c3->Divide(4,4);  
  for( int i=0; i<chRPCR.size(); i++){
    c3->cd(i+1);
    htotR[i]->Draw("");
    htotR[i]->GetYaxis()->SetLabelSize(0.05);
    htotR[i]->GetXaxis()->SetLabelSize(0.07);
    htotR[i]->GetXaxis()->SetTitleSize(0.06);
    htotR[i]->GetXaxis()->SetTitleOffset(0.9);
    hltdcR[i]->GetXaxis()->SetNdivisions(505);
    hltdcR[i]->GetYaxis()->SetNdivisions(505);
  }

  TCanvas *c4 = new TCanvas("c4","c4",3,2,1200,900);
  c4->Divide(4,3);

  c4->cd(1);
  hHitChL->Draw("hbar");
  hHitChL->SetFillColor(kBlue+2);
  hHitChL->GetYaxis()->SetLabelSize(0.05);
  hHitChL->GetXaxis()->SetLabelSize(0.07); 
  hHitChL->GetXaxis()->SetTitleSize(0.06);
  hHitChL->GetXaxis()->SetTitleOffset(0.9);  
  hHitChL->GetYaxis()->SetNdivisions(505);
  
  c4->cd(2);
  hHitChR->Draw("hbar");
  hHitChR->SetFillColor(kBlue+2);
  hHitChR->GetYaxis()->SetLabelSize(0.05);
  hHitChR->GetXaxis()->SetLabelSize(0.07); 
  hHitChR->GetXaxis()->SetTitleSize(0.06);
  hHitChR->GetXaxis()->SetTitleOffset(0.9);
  hHitChR->GetYaxis()->SetNdivisions(505);
  /*
  TImage *img = TImage::Create();
  img->FromPad(c0);
  c0->cd(2);
  img->Flip(90);
  img->Draw();
  */

  c4->cd(5);
  hCh_L1vsL2->Draw("colz");
  hCh_L1vsL2->GetYaxis()->SetLabelSize(0.05);
  hCh_L1vsL2->GetXaxis()->SetLabelSize(0.07); 
  hCh_L1vsL2->GetXaxis()->SetTitleSize(0.06);
  hCh_L1vsL2->GetXaxis()->SetTitleOffset(0.9);
  hCh_L1vsL2->SetStats(0);

  c4->cd(6);
  hCh_R1vsR2->Draw("colz");
  hCh_R1vsR2->GetYaxis()->SetLabelSize(0.05);
  hCh_R1vsR2->GetXaxis()->SetLabelSize(0.07); 
  hCh_R1vsR2->GetXaxis()->SetTitleSize(0.06);
  hCh_R1vsR2->GetXaxis()->SetTitleOffset(0.9);
  hCh_R1vsR2->SetStats(0);


  c4->cd(7);
  hCh_LvsR->Draw("colz");
  hCh_LvsR->GetYaxis()->SetLabelSize(0.05);
  hCh_LvsR->GetXaxis()->SetLabelSize(0.07); 
  hCh_LvsR->GetXaxis()->SetTitleSize(0.06);
  hCh_LvsR->GetXaxis()->SetTitleOffset(0.9);
  hCh_LvsR->SetStats(0);
  hCh_LvsR->SetTitle("Correlation of L and R");

  c4->cd(8);
  hMltptyL->Draw("box");
  hMltptyR->Draw("same");
  hMltptyL->SetStats(0);
  hMltptyL->SetLineColor(kBlue);
  hMltptyR->SetLineColor(kCyan+2);
  hMltptyL->GetXaxis()->SetNdivisions(5);
  hMltptyL->GetYaxis()->SetNdivisions(505);
  TLegend *lg = new TLegend(0.7,0.7,0.9,0.9);
  lg->AddEntry(hMltptyL, "Left", "l");
  lg->AddEntry(hMltptyR, "Right", "l");
  lg->Draw();
  hMltptyL->SetTitle("Hist of Multiplicity");
  
  for (int i=0; i<nTrig; i++){
    c4->cd(9+i);
    hltdc_trig[i]->Draw("box");
    hltdc_trig[i]->GetYaxis()->SetLabelSize(0.05);
    hltdc_trig[i]->GetXaxis()->SetLabelSize(0.07);
    hltdc_trig[i]->GetXaxis()->SetTitleSize(0.06);
    hltdc_trig[i]->GetXaxis()->SetTitleOffset(0.9);
    hltdc_trig[i]->GetYaxis()->SetNdivisions(505);
    hltdc_trig[i]->GetXaxis()->SetNdivisions(505);
  }

  cout << "trigger event = " << trev << endl;
  cout << "event num = " << ev << endl;
  cout << "left event num = " << evL << endl;
  cout << "right event num = " << evR << endl;
  cout << "single strip efficiency L = " << sevL << "/" << trev << " = " << (double) sevL/trev*100 << " ％" << endl;
  cout << "efficiencyL = " << evL << "/" << trev << " = " << (double) evL/trev*100 << "％" << endl;
  cout << "single strip efficiency R = "<< sevR << "/" << trev << " = " << (double) sevR/trev*100 << " ％" << endl;
  cout << "efficiencyR = " << evR << "/" << trev << " = " << (double) evR/trev*100 << "％" << endl;
  cout << "efficiency (several strips) = " << ev << "/" << trev << " = " << (double) ev/trev*100 << "％" << endl;
}

void setStyle() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t Red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t Green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t Blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, Red, Green, Blue, NCont);
    gStyle->SetNumberContours(NCont);
    //  gStyle->SetOptStat(0);
    //gStyle->SetPadGridX(true);
    //gStyle->SetPadGridY(true);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132,"XYZ");
    gStyle->SetLabelFont(132,"XYZ");
    gStyle->SetStatH(6);
    gStyle->SetStatW(0.30);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetLabelSize(0.05,"XY");
    gStyle->SetTitleSize(0.05,"XY");
    gStyle->SetTitleOffset(0.7,"X");
    gStyle->SetTitleOffset(1.3,"X");
    gStyle->SetTitleSize(0.05);
}

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
}