// This code reads a .root input file and extracts the event information. A clustering algorithm is run on the entries.
// The clustering algorithm works as follows:
// 1. For each particle event, pmt hits are sorted in descending order by charge.
// 2. The entry with the highest charge is taken, and pmt hits in the vicinity are combined into a single cluster.
// 3. The algorithm is repeated until all clusters are found.
// 4. The x, y and z positions of the first 10 clusters are plotted, with different colours representing different clusters.
// 5. Cluster information is written on the terminal - eg. total charge in cluster, timing, number of pmt hits in a cluster
//
// The following cuts are applied: Clusters must have more than 3 pmt hits. All pmt hits taken to form a cluster must have at least 5% of the charge of the highest charge entry
// in the cluster. The 'vicinity' is taken as a cuboid with dimensions 355.1m * 355.1m * 559.7m.
//
// The code is run as root -l -x 'ODClusterFeb22.C("&&&&&")' where &&&&& refers to the location of the input file.
//
// Currently, the input file I am using is ../1TeVUpMuJobs_3inchPMTsOD_042percentCoverage_0.root
//
// The radius of HyperK is 3551m and the height is 5597m

#include <stdio.h>     
#include <stdlib.h>   
#include <vector>


// ****************  Opening input file  **************************



void ODClusterFeb22(char *filename=NULL) {

  gROOT->Reset();
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.rootmap");
  }else{
    std::cout << "Can't load WCSim ROOT dictionaries" << std::endl;
  }
  gStyle->SetOptStat(1);

  TFile *f;
  char fTest[128];
  if (filename==NULL){
    std::cout << "Please provide filename in option" << std::endl;
    std::cout << "Will load auto wcsim.root in WCSIMDIR ..." << std::endl;
    char* name = "wcsim.root";
    strncpy(fTest, wcsimdirenv, sizeof(fTest));
    strncat(fTest, "/", (sizeof(fTest) - strlen(fTest)) );
    strncat(fTest, name, (sizeof(fTest) - strlen(fTest)) );
    f = new TFile(fTest);
  }else{
    f = new TFile(filename);
  }
  if (!f->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }else{
    if (filename==NULL) cout << "File open bro: " << fTest << endl;
    else cout << "File open bro: " << filename << endl;
  }


 // **************  Loading the tree  ************************************


  TTree  *wcsimT = (TTree*)f->Get("wcsimT");

  WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();
  wcsimT->SetBranchAddress("wcsimrootevent_OD",&wcsimrootsuperevent);

  // Force deletion to prevent memory leak when issuing multiple
  // calls to GetEvent()
  wcsimT->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);

  const long unsigned int nEntries = wcsimT->GetEntries();
  cout << "Nb of entries " << wcsimT->GetEntries() << endl;

  TTree *wcsimGeoT = (TTree*)f->Get("wcsimGeoT");

  WCSimRootGeom *wcsimrootgeom = new WCSimRootGeom();
  wcsimGeoT->SetBranchAddress("wcsimrootgeom",&wcsimrootgeom);
  wcsimGeoT->GetEntry(0); //Geometry is the same for all events, so can just use GetEntry(0)

  int event;
  float vtxX, vtxY, vtxZ; // vertex coordinate
  float dirX, dirY, dirZ; // particle momentum direction
  float phi, theta; // angles of the particle's momentum
  float energy;

  /*  not sure if i need this   TTree *outTree = new TTree("simulation", "simulation");

  // Set the branches to be saved 
  outTree->Branch("event", &event, "event/I");
  outTree->Branch("vtxX", &vtxX, "vtxX/F");
  outTree->Branch("vtxY", &vtxY, "vtxY/F");
  outTree->Branch("vtxZ", &vtxZ, "vtxZ/F");
  outTree->Branch("dirX", &dirX, "dirX/F");
  outTree->Branch("dirY", &dirY, "dirY/F");
  outTree->Branch("dirZ", &dirZ, "dirZ/F");
  outTree->Branch("phi", &phi, "phi/F");
  outTree->Branch("theta", &theta, "theta/F");
  outTree->Branch("energy", &energy, "energy/F");
 */


  std::vector< vector< TNtuple*>> ClusterContainer; 
  // This is the main container for which all clusters in all events are stored 
  std::vector <TNtuple*> PrimaryParticleInfo;
  // This is the container for primary particle track information
  const float radiusHyperK = 3551;
  const float heightHyperK = 5597;


  // *********  Extracting information from TTree event by event **********
 
  for(long unsigned int iEntry = 0; iEntry <=10; iEntry++){ // I capped it to 10 instead of nEntries for now, just for ease of plotting.)
      wcsimT->GetEvent(iEntry); 
      std::cout<< "Reading event " << iEntry <<std::endl;
      // Point to event iEntry inside WCSimTree

      // Initialising variables for each event loop
      const long unsigned int nTriggers = wcsimrootsuperevent->GetNumberOfEvents();
      const long unsigned int nSubTriggers = wcsimrootsuperevent->GetNumberOfSubEvents();
      TNtuple *ntuple = new TNtuple("ntuple","x y z charge ntuple","pmtX:pmtY:pmtZ:cylLoc:charge:time"); 
      TNtuple *trackinfo = new TNtuple("trackinfo","trackinfo","vtxX:vtxY:vtxZ:dirX:dirY:dirZ:energy:theta:phi"); // might need to set the branch address at some point
 
      Double_t maxCharge = 0.0;
      Double_t maxpmtX = 0.0;
      Double_t maxpmtY = 0.0;
      Double_t maxpmtZ = 0.0;
      Double_t maxtime = 0.0;
      int maxcylLoc = 0;
      Double_t nClusters = 0;
      std::vector <TNtuple*> SingleEventClusters;

      for(long unsigned int iTrig = 0; iTrig < nTriggers; iTrig++){
	 WCSimRootTrigger *wcsimroottrigger = wcsimrootsuperevent->GetTrigger(iTrig);
	 int ncherenkovdigihits = wcsimroottrigger->GetNcherenkovdigihits();

      
  // ******* extracting primary particle information ***************

	 if (iTrig == 0){
	   int numTracks = wcsimroottrigger->GetNtrack();
	   WCSimRootTrack *track = (WCSimRootTrack* ) wcsimroottrigger->GetTracks()->At(0);
	   
	   vtxX = wcsimroottrigger->GetVtx(0);
	   vtxY = wcsimroottrigger->GetVtx(1);
	   vtxZ = wcsimroottrigger->GetVtx(2);
	   dirX = track->GetDir(0);
	   dirY = track->GetDir(1);
	   dirZ = track->GetDir(2);
	   energy = track->GetE();
	   theta = abs( atan (dirX/dirZ));
	   phi = abs( atan( dirY/dirX) );
	   trackinfo->Fill(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,energy,theta,phi);
	   //trackinfo-> Scan();

	   float vX, vY, vZ, dX,dY, dZ, Energy, Theta, Phi;
	   trackinfo->SetBranchAddress("vtxX",&vX);
	   trackinfo->SetBranchAddress("vtxY",&vY);
	   trackinfo->SetBranchAddress("vtxZ",&vZ);
	   trackinfo->SetBranchAddress("dirX",&dX);
	   trackinfo->SetBranchAddress("dirY",&dY);
	   trackinfo->SetBranchAddress("dirZ",&dZ);
	   trackinfo->SetBranchAddress("energy",&Energy);
	   trackinfo->SetBranchAddress("theta",&Theta);
	   trackinfo->SetBranchAddress("phi",&Phi);
	   

	   
	   PrimaryParticleInfo.push_back(trackinfo);
	 }


  // ****** Filling up a container with all the pmt hits from a single event ******

	 for (int i=0; i<ncherenkovdigihits; i++){
	   WCSimRootCherenkovDigiHit *hit = (WCSimRootCherenkovDigiHit*)
	     (wcsimroottrigger->GetCherenkovDigiHits()->At(i));
	   
	   int tubeID = hit->GetTubeId()-1 ; //not sure why -1. need to ask stephane. took it from his email.
	   WCSimRootPMT pmt = wcsimrootgeom->GetPMT(tubeID);
	   Double_t cylLoc = pmt.GetCylLoc();
	   Double_t pmtX = pmt.GetPosition(0);
	   Double_t pmtY = pmt.GetPosition(1);
	   Double_t pmtZ = pmt.GetPosition(2);
	   Double_t charge = hit->GetQ();
	   Double_t time = hit->GetT();
	   ntuple->Fill(pmtX,pmtY,pmtZ,cylLoc,charge,time);
	  
	   // Finding pmt hit with maximum charge in the event
	   if (charge>maxCharge){
	     maxcylLoc = cylLoc;
	     maxCharge = charge;
	     maxpmtX = pmtX;
	     maxpmtY = pmtY;
	     maxpmtZ = pmtZ;
	     maxtime = time;
	   }
	 }
      }


  // ************* Clustering algorithm *****************
       
      do {     
 
	float  pX,pY,pZ, Q, T, cL;
	TNtuple *maxchargentuple = new TNtuple("maxchargentuple","x y z charge ntuple","pmtX:pmtY:pmtZ:cylLoc:charge:time");
	TNtuple *n2ntuple = new TNtuple("n2ntuple","x y z charge ntuple","pmtX:pmtY:pmtZ:cylLoc:charge:time");
	ntuple->SetBranchAddress("pmtX",&pX);
	ntuple->SetBranchAddress("pmtY",&pY);
	ntuple->SetBranchAddress("pmtZ",&pZ);
	ntuple->SetBranchAddress("cylLoc",&cL);
	ntuple->SetBranchAddress("charge",&Q);
	ntuple->SetBranchAddress("time",&T);
	n2ntuple->SetBranchAddress("pmtX",&pX);
	n2ntuple->SetBranchAddress("pmtY",&pY);
	n2ntuple->SetBranchAddress("pmtZ",&pZ);
	n2ntuple->SetBranchAddress("cylLoc",&cL);
	n2ntuple->SetBranchAddress("charge",&Q);
	n2ntuple->SetBranchAddress("time",&T);
	maxchargentuple->SetBranchAddress("pmtX",&pX);
	maxchargentuple->SetBranchAddress("pmtY",&pY);
	maxchargentuple->SetBranchAddress("pmtZ",&pZ);
	maxchargentuple->SetBranchAddress("cylLoc",&cL);
	maxchargentuple->SetBranchAddress("charge",&Q);
	maxchargentuple->SetBranchAddress("time",&T);
   
	Double_t maxnextCharge = 0.0;
	Double_t maxnextpmtX = 0.0;
	Double_t maxnextpmtY = 0.0;
	Double_t maxnextpmtZ = 0.0;
	Double_t  maxnextcylLoc = 0;
	Double_t maxnexttime = 0.0;
	Double_t totalqincluster = 0.0; 

	for ( Long64_t ient = 0; ient < ntuple -> GetEntries(); ient++ ) {
        
	  ntuple->GetEntry(ient);
	  
	  if ((pX>maxpmtX-(0.2*radiusHyperK)) && (pX<maxpmtX+(0.2*radiusHyperK)) && (pY>maxpmtY-(0.2*radiusHyperK)) && (pY<maxpmtY+(0.2*radiusHyperK)) && (pZ>maxpmtZ-(0.1*heightHyperK)) && (pZ<maxpmtZ+(0.1*heightHyperK) ) && (Q>(0.05*maxCharge) ) ){
	    //using 10% diameter (2*10% radius) and 10% height

	    maxchargentuple->Fill(pX,pY,pZ,cL,Q,T);
	    totalqincluster = totalqincluster+Q;
	  }
	  
     
	  else if ((pX>maxpmtX-(0.2*radiusHyperK)) && (pX<maxpmtX+(0.2*radiusHyperK)) && (pY>maxpmtY-(0.2*radiusHyperK)) && (pY<maxpmtY+(0.2*radiusHyperK)) && (pZ>maxpmtZ-(0.1*heightHyperK)) && (pZ<maxpmtZ+(0.1*heightHyperK) ) && (Q<(0.05*maxCharge) ) ){                       
	    continue;       
	  }
	  
	  else { n2ntuple->Fill(pX,pY,pZ,cL,Q,T); 

	    if (Q>maxnextCharge){                        // This is the highest charge outside of the first cluster, and its associated coordinates
	      maxnextCharge = Q;
	      maxnextpmtX = pX;
	      maxnextpmtY = pY;
	      maxnextpmtZ = pZ;
	      maxnextcylLoc = cL;
	      maxnexttime = T;
	    }
	  }
	  
	}
	
	SingleEventClusters.push_back(maxchargentuple);
	nClusters++;
	cout << "Max charge event in cluster " << nClusters << " is: " << " X: " << maxpmtX << " Y: " <<maxpmtY << " Z: " <<maxpmtZ << " Q: " <<  maxCharge << " time: " <<  maxtime << endl;
	Double_t nHitsInEachCluster = maxchargentuple->GetEntries();
	cout << "Total charge in cluster " << nClusters << " is " << totalqincluster << endl;
	cout << "Number of PMT hits in cluster " << nClusters << "  is: " << nHitsInEachCluster << endl;

	TNtuple *ntuple = new TNtuple("ntuple","x y z charge ntuple","pmtX:pmtY:pmtZ:cylLoc:charge:time"); // deleting previous ntuple and using new, reduced ntuple for next clustering loop.
	ntuple = n2ntuple; 
	maxCharge = maxnextCharge;
	maxpmtX = maxnextpmtX;
	maxpmtY = maxnextpmtY;
	maxpmtZ = maxnextpmtZ;
	maxtime=maxnexttime;
	maxcylLoc=maxnextcylLoc;
      }
      while (nHitsInEachCluster > 3 &&  iEntry <= nEntries );

      cout << "Event " << iEntry << " has " << nClusters << " clusters." << endl;
      ClusterContainer.push_back(SingleEventClusters);

   }


   // ************* Drawing *********************************


   TCanvas *c1;
   c1 = new TCanvas("ntuple","ntuple", 2400,1200);
   
   c1->Divide(5,2);

   for (long unsigned int nEvents = 0; nEvents < ClusterContainer.size(); nEvents++){

     c1->cd(nEvents+1);      
     TGeoTube *myCy1 = new TGeoTube(0,3551,2798.5);        
     myCy1->Draw();


     // ******* Drawing particle track ******************** //
     PrimaryParticleInfo[nEvents]->GetEntry(0);
     cout << "vtxX " << vX << " vtxY " << vY << endl;
     TPolyLine3D *pl3d1 = new TPolyLine3D(2);
   // set points - I manually set this for now, but want to find a more efficient way to set this from data.
     pl3d1->SetPoint(0, vX,vY,vZ);
     pl3d1->SetPoint(1, vX,vY,3500);
     pl3d1->SetLineWidth(3);
     pl3d1->SetLineColor(2);
     pl3d1->Draw("same");


     // ********* Drawing pmt hits ***********************     
     for (long unsigned int eachCluster = 0; eachCluster < ClusterContainer[nEvents].size(); eachCluster++){

       ClusterContainer[nEvents][eachCluster]->SetMarkerStyle(21);
       ClusterContainer[nEvents][eachCluster]->SetMarkerColor(eachCluster+1);  // red
       ClusterContainer[nEvents][eachCluster]->Draw("pmtZ:pmtY:pmtX","", "same");

       
     }
     
   }
   

   /*
   TH1F *h1 = new TH1F("h1","hist from tree",50, -3300, 3300);
   ClusterContainer[1][1]->Draw("charge:pmtX>>h1"); 
   h1->Fit("gaus","same");
   */
}
