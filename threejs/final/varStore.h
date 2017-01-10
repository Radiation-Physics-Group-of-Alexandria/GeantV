class varStore{
	public:
		bool MapVol[100000];
		int kVolume=
		0;
		int nVolumes=0;
		int clonedVolumes=0;
		int MaxVisiLevel=999;
		bool OneCopy=true;
		int level=0;
		int previousLevel;
		int kgeom=-1;
		int nColors=0;
		int nsteps=0;
		int nTotalNodes=0;
		TGeoVolume *top; 
		TGeoNode *current;
		TGeoShape *shape;
		TGeoMatrix *matrix;
		Double_t Mtrx[16];
		Double_t nM=0;
		TString Comp="";
		ofstream PrimFile;
		ofstream ClonesFile;
		TString PrimFilename="primitives.js";
		TString ClonesFilename="completetree.js";
		int nMaxVol=0;

		int narcseg=20;
    	int ntubseg=10;
    	int ncseg=45; //number of segments on the circonference
    	int nzseg=45; //number of segments on Z
    	int nphiseg=45;
    	int ntheseg=45;
    	int geoId=-1;
};