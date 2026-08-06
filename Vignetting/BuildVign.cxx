//Use ./prog <rootfile>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include <algorithm>

using namespace std;

void findmax(shared_ptr<TH2D>,int& ,int&);

int main(int argc, char** argv)
{
    string namefile= argv[1];
    string namepic= "pic_0";
    constexpr int bins = 3001;
    constexpr int binx = 4096;
    constexpr int biny = 2304;
    constexpr int centreQinx = 2025;                    //Quest Xenon lens (david) 2025     Quest EHD lens ap=0.85 (stefano) 1968   Quest EHD lens ap=0.95 (stefano) 1941
    constexpr int centreQiny = 914;                    //Quest Xenon lens (david) 914     Quest EHD lens ap=0.85 (stefano) 1019    Quest EHD lens ap=0.95 (stefano) 1000
    constexpr int centreQoutx = 2048;
    constexpr int centreQouty = 1152;

    auto h1vign = make_shared<TH1D>("Vign1D","Vign1D",bins,-0.5,3000.5);
    auto normmap = make_shared<TH2F>("normmap","normmap",binx,0,binx,biny,0,biny);
    auto hcont = make_shared<TH1D>("hcon","hcon",bins,-0.5,3000.5);
    auto fin = shared_ptr<TFile>(TFile::Open(namefile.c_str(),"READ"));
    auto h2vign = shared_ptr<TH2D>((TH2D*) fin->Get(namepic.c_str()));
    auto fout = shared_ptr<TFile>(TFile::Open("VignQ1D.root","RECREATE"));

    //Find average max value in the centre of the input image
    double max=0;
    for(int i=centreQinx-3;i<=centreQinx+3;i++)
    {
        for(int j=centreQiny-3;j<=centreQiny+3;j++)
        {
            max+=h2vign->GetBinContent(i,j);
        }
    }
    max/=49;
    //Fill the 1D histogram with the relative intensity of the pixel at a distance r from centre (if more pixels are at a distance r they will be summed and the average will be taken (next loop)
    for(int i=1;i<=binx;i++)
    {
        for(int j=1;j<=biny;j++)
        {
            double r=sqrt((i-centreQinx)*(i-centreQinx)+(j-centreQiny)*(j-centreQiny));
            h1vign->Fill(r,h2vign->GetBinContent(i,j)/max);
            hcont->Fill(r);
        }
    }
    for(int i=1;i<=bins;i++)
    {
        if(hcont->GetBinContent(i)!=0)
        {
            h1vign->SetBinContent(i,h1vign->GetBinContent(i)/hcont->GetBinContent(i));
            h1vign->SetBinError(i,sqrt(h1vign->GetBinContent(i)*hcont->GetBinContent(i))/hcont->GetBinContent(i));      //Technically this makes no sense. I just need this to weight the vignetting values for the number of entries
        }
        
    }
    //Fill the 2D map of vignetting calculating the distance of the pixels from the new centre
    for(int i=1;i<=binx;i++)
    {
        for(int j=1;j<=biny;j++)
        {
            double r=sqrt((i-centreQoutx)*(i-centreQoutx)+(j-centreQouty)*(j-centreQouty));
            int rbin= int(r);
            if(r-rbin>=0.5) rbin++;
            rbin++;     //this because in ROOT histogram used bin 1 has  edges from -0.5 to 0.5
            //if(i%1000==0) {if(j%1000==0) cout<< r << " "<<h1vign->GetBinContent(r+1) <<" "<< h1vign->GetBinContent(1349) <<endl;}
            normmap->SetBinContent(i,j,h1vign->GetBinContent(rbin));
        }
    }
    //Save in the fullinfo rootfile
    h1vign->Rebin(4);
    h1vign->Scale(1./h1vign->GetBinContent(1));
    normmap->SetName("normmap_full");

    double max2=0;
    int centercount=0;
    for(int i=centreQoutx-3;i<=centreQoutx+3;i++)                                   //49 pixels for QUEST if not rebinned (xc-3, xc+3)      9 for 4x4 rebin (equivalent to 144 original pixels)
    {
        for(int j=centreQouty-3;j<=centreQouty+3;j++)
        {
            max2+=normmap->GetBinContent(i,j);
            centercount++;
        }
    }
    max2/=centercount;
    cout<<normmap->GetBinContent(normmap->GetMaximumBin())<<" "<<max2<<endl;
    normmap->Scale(1./max2);
    normmap->Write();
    h1vign->Write();

    //Rebin and normalise and store in the root file for the reco
    auto foutreco = shared_ptr<TFile>(TFile::Open("VignQ1D_reco.root","RECREATE"));
    normmap->SetName("normmap");
    normmap->Rebin2D(4,4);
    //Findmax in new vignett
    max2=0;
    centercount=0;
    constexpr int centreQoutx_rebin=centreQoutx/4;
    constexpr int centreQouty_rebin=centreQouty/4;
    for(int i=centreQoutx_rebin-1;i<=centreQoutx_rebin+1;i++)                                   //49 pixels for QUEST if not rebinned (xc-3, xc+3)      9 for 4x4 rebin (equivalent to 144 original pixels)
    {
        for(int j=centreQouty_rebin-1;j<=centreQouty_rebin+1;j++)
        {
            max2+=normmap->GetBinContent(i,j);
            centercount++;
        }
    }
    max2/=centercount;
    normmap->Scale(1./max2);
    //Adjust close to 1 pixels to 1
    for(int i=centreQoutx_rebin-45;i<centreQoutx_rebin+45;i++)                                //401x401 pixels for no rebin           91x91 pixels for rebin 
    {
        for(int j=centreQouty_rebin-45;j<centreQouty_rebin+45;j++)
        {
            double r=sqrt((i-centreQoutx_rebin)*(i-centreQoutx_rebin)+(j-centreQouty_rebin)*(j-centreQouty_rebin));
            if(normmap->GetBinContent(i,j)>0.99 && r<25) normmap->SetBinContent(i,j,1);            //r<100 for no rebin        25 for rebin4
        }
    }

    normmap->Write();
    return 0;
}

void findmax(shared_ptr<TH2D> h2,int& x,int& y)
{
    double maxval=0;
    for(int i=1;i<4097;i++)
    {
        for(int j=1;j<2305;j++)
        {
            double te=h2->GetBinContent(i,j);
            if(te>maxval)
            {
                maxval=te; x=i; y=j;
            }
        }
    }
    return;
}