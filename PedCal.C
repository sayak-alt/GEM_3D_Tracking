#include "TFile.h"
#include "TH1F.h"
#include <assert.h>
#include <utility>
#include "TF1.h"
#include "TChain.h"
#include "TMath.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <map>

#include <TCanvas.h>

Int_t create_pedDB=0;
Double_t fitf(Double_t *v, Double_t *par)
{
  	Double_t arg = 0;
  	if (par[2] != 0) arg = (v[0] - par[1])/par[2];
	Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
   	return fitval;
}


void PedCal()
{
	Int_t NoOfDet=1;
	ofstream fh_out2;
    //run number for the pedestal file
	Int_t argc=3910;
	
    const int NDet = 1.0;

    	vector<int> xstrip;
    	vector<int> ystrip;
    	int dumy1,dumy2;
    	for(int pq=0;pq<NoOfDet;pq++)
    	{
            //change the no of strips are there in the output depending on how many APV's are being used
            dumy1 = 256; // X-axis
            dumy2 = 512-128; //Y-axis GEM-I
         //  dumy2 = 512-128; //Y-axis GEM-II & GEM-III
		cout<<"xch "<<dumy1<<" ych "<<dumy2<<endl;
    		xstrip.push_back(dumy1);
    		ystrip.push_back(dumy2);
    	}
    

	cout<<"Run number "<<argc<<endl;
	TChain *T = new TChain("T");
    //Pedestal root file path
    T->Add(Form("/Users/sayakchatterjee/Desktop/Sam_ug_2025/data_files_root/test_%d.root",argc));
  //  T->Add(Form("/Users/sayakchatterjee/Desktop/Sam_ug_2025/data_files_root/pedestal.root",argc));
    

	Double_t* x_strip=new Double_t[NDet];
        Double_t* y_strip=new Double_t[NDet];//# of strips in X/Y
        Double_t*** xadc = new Double_t**[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
        {
                //const int xch = (int) xstrip[detnumber];
                xadc[detnumber]=new Double_t*[6];              
                for(int NoOfSample=0;NoOfSample<6;NoOfSample++)
                        xadc[detnumber][NoOfSample]=new Double_t[xstrip[detnumber]];
        }

        Double_t*** yadc = new Double_t**[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
        {
                yadc[detnumber]=new Double_t*[6];
                for(int NoOfSample=0;NoOfSample<6;NoOfSample++)
                        yadc[detnumber][NoOfSample]=new Double_t[ystrip[detnumber]];
        }

        Double_t** xstripID = new Double_t*[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
                xstripID[detnumber]=new Double_t[xstrip[detnumber]];

        Double_t** ystripID = new Double_t*[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
                ystripID[detnumber]=new Double_t[ystrip[detnumber]];

        Double_t** xped= new Double_t*[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
                xped[detnumber]=new Double_t[xstrip[detnumber]];

        Double_t** xrms= new Double_t*[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
                xrms[detnumber]=new Double_t[xstrip[detnumber]];

        Double_t** yped= new Double_t*[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
                yped[detnumber]=new Double_t[ystrip[detnumber]];

        Double_t** yrms= new Double_t*[NDet];
        for(int detnumber=0;detnumber<NDet;detnumber++)
                yrms[detnumber]=new Double_t[ystrip[detnumber]];

	map<int, vector<TH1F*> > hx;
	map<int, vector<TH1F*> > hy;

	for(int detnumber=0;detnumber<NoOfDet;detnumber++)
    	{
//change the "+1" to "+2" for GEM-II and "+3" for GEM-III to access the respective data branch
		T->SetBranchAddress(Form("sbs.gems.x%d.nch",detnumber+1), &x_strip[detnumber]);
    		T->SetBranchAddress(Form("sbs.gems.y%d.nch",detnumber+1), &y_strip[detnumber]);
	    	T->SetBranchAddress(Form("sbs.gems.x%d.strip",detnumber+1), xstripID[detnumber]);
	    	T->SetBranchAddress(Form("sbs.gems.y%d.strip",detnumber+1), ystripID[detnumber]);
		for(int ij=0;ij<6;ij++)
		{
	      		T->SetBranchAddress(Form("sbs.gems.x%d.adc%d",detnumber+1,ij),xadc[detnumber][ij]);
	      		T->SetBranchAddress(Form("sbs.gems.y%d.adc%d",detnumber+1,ij),yadc[detnumber][ij]);
	  	}
		
		for(int ij=0;ij<xstrip[detnumber];ij++)
	    	{
	       		hx[detnumber].push_back(new TH1F(Form("hx_%d_%d",detnumber+3,ij+3),Form("Det %d Pedestal x_%d",detnumber+3,ij+3),3000,0,3000));
	    	}
		
		for(int ij=0;ij<ystrip[detnumber];ij++)
	    	{
	       		hy[detnumber].push_back(new TH1F(Form("hy_%d_%d",detnumber+3,ij+3),Form("Det %d Pedestal y_%d",detnumber+3,ij+3),3000,0,3000));
	    	}
	}
	    
	int entries = T->GetEntries();
	    	
	for(int ij=0;ij<entries;ij++)
	
	{
		T->GetEntry(ij);
		for(int detnumber=0;detnumber<NDet;detnumber++)
		{
			Double_t dummy = 0;  
			for(int pq=0;pq<xstrip[detnumber];pq++)
	    		{
	    			dummy=0;
	    			for(int ik=0;ik<6;ik++)
	    			{
	    				dummy += xadc[detnumber][ik][pq];
	    			}
	    			dummy = dummy/6.;
	    			hx[detnumber][pq]->Fill(dummy);
			}
	    		for(int pq=0;pq<ystrip[detnumber];pq++)
	    		{
	    			dummy = 0;
	    			for(int ik = 0;  ik<6; ik++)
	    			{
	    				dummy += yadc[detnumber][ik][pq];
	    			}
	    
				dummy = dummy/6;
	    			hy[detnumber][pq]->Fill(dummy);
	    		}
	    	}
	}
	    	
    	ofstream fh_out;
        //Change the output file names. We need different output files for different GEMs
        TString outfile="SelectPulse_3910_x1y1.txt";//for GEM-I
 //       TString outfile="SelectPulse_3417_x2y2.txt";//for GEM-II
 //   TString outfile="SelectPulse_3417_x3y3.txt";//for GEM-III


    	fh_out.open(outfile);
    	if(!fh_out)
    	{
    		cout<<"Cannot open output file " <<outfile<<endl;
    		return;
    	}

    	ofstream fh_out1;
        TString outfile1="SelectPulse_ped_x3y3_additional.txt"; //for GEM-I
//          TString outfile1="SelectPulse_3417_x2y2_additional.txt"; //for GEM-II
//              TString outfile1="SelectPulse_3417_x3y3_additional.txt"; //for GEM-III


    	fh_out1.open(outfile1);
    	if(!fh_out1)
    	{
    		cout<<"Cannot open output file " <<outfile1<<endl;
    		return;
    	}
        if(create_pedDB==1){	
    	TString dbfile="db_sbs.gems.dat";
    	fh_out2.open(dbfile,std::ofstream::out|std::ofstream::app);
    	if(!fh_out2)
    	{
    		cout<<"Cannot open output file " <<dbfile<<endl;
    		return;
    	}
	}
	Double_t xcentroid;
	Double_t ycentroid;
	Double_t xsigma;
	Double_t ysigma;
	Double_t xamp;
	Double_t yamp;

	map<int, vector<TCanvas*> > plotx;
    	map<int, vector<TCanvas*> > ploty;
	for(int detnumber=0;detnumber<NDet;detnumber++)
	{
    		int xcanvasIndex=0;
		for(int ij=0; ij<xstrip[detnumber];ij++)
    		{
    			hx[detnumber][ij]->GetXaxis()->SetRangeUser(1,3000);
	    		xcentroid=(hx[detnumber][ij]->GetMaximumBin()-1);
    			xamp=(hx[detnumber][ij]->GetBinContent(xcentroid));
    			xsigma=30;
//                xsigma=1; //testing new sigma values
	    		TF1 *func = new TF1("fit",fitf,xcentroid-200,xcentroid+200,3);
    			/*if(ij%125==0)
			cout<<"  xch  "<<ij<<" xcentroid  "<<xcentroid[ij]<<
	    		"  xamp "<<xamp[ij]<<"  xsigma  "<<xsigma[ij]<<endl;*/
    			func->SetParameters(xamp,xcentroid,xsigma);
    			func->SetParNames("Constant","Mean_value","Sigma");
	    		hx[detnumber][ij]->Fit("fit","RQ0","",xcentroid-200,xcentroid+200);
    
			fh_out<<fixed<<setprecision(2)<<func->GetParameter(1)<<"\t"<<fixed<<setprecision(2)<<fabs(func->GetParameter(2))*sqrt(6)<<endl;
			xped[detnumber][ij]=func->GetParameter(1);
			xrms[detnumber][ij]=fabs(func->GetParameter(2));
    			if(ij%125==0)
			{
				cout<<"After fitting for Det "<<detnumber<<"  xcentroid "<<func->GetParameter(1)<<"  xamplitude  "<<func->GetParameter(0)<<" xsigma  "<< func->GetParameter(2)*sqrt(6)<<endl;
	    			plotx[detnumber].push_back(new TCanvas(Form("xplot_%d_%d",detnumber+3,ij+3),Form("Det %d Ped x strip %d",detnumber+3,ij+3),800,800));
	    			plotx[detnumber][xcanvasIndex]->SetLogy();
		    		hx[detnumber][ij]->GetXaxis()->SetTitleSize(0.04);
		    		hx[detnumber][ij]->GetXaxis()->SetTitleOffset(1.0);
	    			hx[detnumber][ij]->GetXaxis()->SetTitle("ADC");
	    			hx[detnumber][ij]->GetYaxis()->SetTitleOffset(1.0);
		    		hx[detnumber][ij]->GetYaxis()->SetTitleSize(0.04);
		    		hx[detnumber][ij]->GetYaxis()->SetTitle("Counts");
				hx[detnumber][ij]->Draw();
	    			func->Draw("same");
		    		xcanvasIndex++;
    			}
    		}
	    	int ycanvasIndex=0;
    		for(int ij=0; ij<ystrip[detnumber];ij++)
	    	{
    			hy[detnumber][ij]->GetXaxis()->SetRangeUser(1,3000);
    			ycentroid=hy[detnumber][ij]->GetMaximumBin()-1;
	    		yamp=hy[detnumber][ij]->GetBinContent(ycentroid);
    			ysigma=30;
//                ysigma=1;//testing new sigma values
    			TF1 *func1 = new TF1("fit1",fitf,ycentroid-200,ycentroid+200,3);
	    		/*if(ij%50==0)
    			cout<<"  ych  "<<ij<<" ycentroid  "<<ycentroid[ij]<<
	        	"  yamp "<<yamp[ij]<<"  ysigma  "<<ysigma[ij]<<endl;*/
	    		func1->SetParameters(yamp,ycentroid,ysigma);
    			func1->SetParNames("Constant","Mean_value","Sigma");
			hy[detnumber][ij]->Fit("fit1","RQ0","",ycentroid-200,ycentroid+200);
    			//cout<<"After fitting  centroid "<<func1->GetParameter(1)<<"  yamplitude  "<<func1->GetParameter(0)<<" ysigma  "<< func1->GetParameter(2)<<endl;
    			fh_out<<fixed<<setprecision(2)<<func1->GetParameter(1)<<"\t"<<fixed<<setprecision(2)<<fabs(func1->GetParameter(2))*sqrt(6)<<endl;
			yped[detnumber][ij]=func1->GetParameter(1);
			yrms[detnumber][ij]=fabs(func1->GetParameter(2));
			//fh_out1<<xstripID[detnumber][ij]<<"\t"<<func1->GetParameter(1)<<"\t"<<fabs(func1->GetParameter(2))<<endl;
		    	if(ij%125==0)
			{
	    			cout<<"After fitting for Det "<<detnumber<<" centroid "<<func1->GetParameter(1)<<"  yamplitude  "<<func1->GetParameter(0)<<" ysigma  "<< func1->GetParameter(2)<<endl;
	    			ploty[detnumber].push_back(new TCanvas(Form("Det %d yplot_%d",detnumber+3,ij+3),Form("Det %d Ped y strip %d",detnumber+3,ij+3),800,800));
		    		ploty[detnumber][ycanvasIndex]->SetLogy();
		    		hy[detnumber][ij]->GetXaxis()->SetTitleSize(0.04);
	    			hy[detnumber][ij]->GetXaxis()->SetTitleOffset(1.0);
	    			hy[detnumber][ij]->GetXaxis()->SetTitle("ADC");
		    		hy[detnumber][ij]->GetYaxis()->SetTitleOffset(1.0);
		    		hy[detnumber][ij]->GetYaxis()->SetTitleSize(0.04);
	    			hy[detnumber][ij]->GetYaxis()->SetTitle("Counts");
	   			hy[detnumber][ij]->Draw();
		    		func1->Draw("same");
	    			ycanvasIndex++;
	    		}
    		}
	}
	fh_out.close();
	
	for(int detnumber=0;detnumber<NDet;detnumber++)
	{
		fh_out1<<Form("\nsbs.gems.x%d.ped = ",detnumber+3);
		if(create_pedDB==1)fh_out2<<Form("\nsbs.gems.x%d.ped = ",detnumber+3);
		for(int ij=0; ij<xstrip[detnumber];ij++)
    		{
			if( ij%8==0)
			{
				fh_out1<<"\\"<<endl;
				if(create_pedDB==1)fh_out2<<"\\"<<endl;
			}
			fh_out1<<fixed<<setprecision(0)<<xstripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<xped[detnumber][ij]<<" ";
			if(create_pedDB==1)fh_out2<<fixed<<setprecision(0)<<xstripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<xped[detnumber][ij]<<" ";
		}
		fh_out1<<Form("\nsbs.gems.x%d.rms = ",detnumber+3);
		if(create_pedDB==1)fh_out2<<Form("\nsbs.gems.x%d.rms = ",detnumber+3);
		for(int ij=0; ij<xstrip[detnumber];ij++)
    		{
			if(ij%8==0)
			{
				fh_out1<<"\\"<<endl;
				if(create_pedDB==1)fh_out2<<"\\"<<endl;
			}
			fh_out1<<fixed<<setprecision(0)<<xstripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<xrms[detnumber][ij]*sqrt(6)<<" ";
			if(create_pedDB==1)fh_out2<<fixed<<setprecision(0)<<xstripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<xrms[detnumber][ij]*sqrt(6)<<" ";
		}
		fh_out1<<Form("\nsbs.gems.y%d.ped = ",detnumber+3);
		if(create_pedDB==1)fh_out2<<Form("\nsbs.gems.y%d.ped = ",detnumber+3);
		for(int ij=0; ij<ystrip[detnumber];ij++)
    		{
			if( ij%8==0)
			{
				fh_out1<<"\\"<<endl;
				if(create_pedDB==1)fh_out2<<"\\"<<endl;
			}
			fh_out1<<fixed<<setprecision(0)<<ystripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<yped[detnumber][ij]<<" ";
			if(create_pedDB==1)fh_out2<<fixed<<setprecision(0)<<ystripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<yped[detnumber][ij]<<" ";
		}
		fh_out1<<Form("\nsbs.gems.y%d.rms = ",detnumber+3);
		if(create_pedDB==1)fh_out2<<Form("\nsbs.gems.y%d.rms = ",detnumber+3);
		for(int ij=0; ij<ystrip[detnumber];ij++)
    		{
			if( ij%8==0)
			{
				fh_out1<<"\\"<<endl;
				if(create_pedDB==1)fh_out2<<"\\"<<endl;
			}
			fh_out1<<fixed<<setprecision(0)<<ystripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<yrms[detnumber][ij]/sqrt(6)<<" ";
			if(create_pedDB==1)fh_out2<<fixed<<setprecision(0)<<ystripID[detnumber][ij]<<" "<<fixed<<setprecision(2)<<yrms[detnumber][ij]*sqrt(6)<<" ";
		}

	}
	fh_out1.close();
	if(create_pedDB==1)fh_out2.close();
	for(int detnumber=0;detnumber<NDet;detnumber++)
	{
		for(int NoOfSample=0;NoOfSample<6;NoOfSample++)
		{
			delete [] xadc[detnumber][NoOfSample];
			delete [] yadc[detnumber][NoOfSample];
		}
		delete [] xadc[detnumber];
		delete [] yadc[detnumber];
		delete [] xped[detnumber];
		delete [] yped[detnumber];
		delete [] xrms[detnumber];
		delete [] yrms[detnumber];
	}
	delete [] xadc;
	delete [] yadc;
}

