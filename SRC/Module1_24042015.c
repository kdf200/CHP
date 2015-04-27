//***********************************************************************
//**               Module 1 of Canopy Height Program                   **
//***********************************************************************
//** Reads PulseWaves data format and converts .pls/ .wvs into         **
//** a .txt file for further use with remaining modules of the program.**
//** Authors: Karolina Fieber, Ian Davenport                           **
//** Adapted from GeoCodeWFExtract program by Ian Davenport            **
//** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **
//** University of Reading, UK              November 2014 - March 2015 **
//***********************************************************************

//#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <malloc.h>


unsigned char readuchar(FILE *f) {
  unsigned char v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

char readchar(FILE *f) {
  char v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

unsigned short int readushort(FILE *f) {
  unsigned short int v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

short int readshort(FILE *f) {
  short int v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

unsigned long readulong(FILE *f) {
  unsigned long v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

long readlong(FILE *f) {
  long v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

unsigned __int64 readi64(FILE *f) {
  unsigned __int64 v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

__int64 readlonglong(FILE *f) {
  __int64 v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

float readfloat(FILE *f) {
  float v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}

double readdouble(FILE *f) {
  double v;
  fread((void*)(&v), sizeof(v), 1, f);
  return v;
}



main()
{
	/* File pointers */

	FILE *wvs_fp, *pls_fp, *param_fp, *out_waveform_fp,*log_fp;

	/* Variable declarations */

	char wvs_filename[200], pls_filename[200], param_filename[200], out_waveform_filename[200], filename[200],log_filename[200]; 
	int area_type;
	double easting_centre, northing_centre, radius;
	double south_edge, north_edge, west_edge, east_edge, low_edge, high_edge;
	int extract_threshold;
	int output_waveform;
	double min_easting, max_easting, min_northing, max_northing, min_elevation, max_elevation;

	///////////////////////////////////////////////////////////////
	//.pls header variables
	char FileSgn[17], SystemID[65], Software[65];
	unsigned char ProjectID4[9], VersionMajor, VersionMinor;
	unsigned long GlobalParam, FileSourceID, ProjectID1, NumberOfVLR, NumberOfAVLR;
	unsigned short ProjectID2, ProjectID3, CreationDayofYear, FileCreationYear, HeaderSize ;
	long long OffsetToPulse, NumberOfPulses, ReservedPLS, MinT, MaxT;
	long PulseFormat, PulseAttributes, PulseSize, PulseCompression;
	double TScaleFactor, TOffset, XScale, YScale, ZScale, XOffset, YOffset, ZOffset, MinX, MaxX, MinY, MaxY, MinZ,  MaxZ;

	// .pls VLR Header variables //200 size is the maximum number of AVLRs/VLRs
	char UserID_VLR[200][17],Description_VLR[200][65];
	unsigned long RecordID_VLR[200],Reserved_VRL[200];
	//unsigned long New_reserved[200];
	long long RecordLenghtAfterHeader_VLR[200];

	// .pls Pulse records
	long long GPSTime, OffsetToWaves;
	long AnchorX, AnchorY, AnchorZ, TargetX, TargetY, TargetZ;
	unsigned short FirstSample, LastSample;
	unsigned char Intensity, Classification,PulseDescrIndex,Reserved_pulse;

	// .pls AVLR header variables //200 size is the maximum number of AVLRs/VLRs
	//char UserID_AVLR[200][17],Description_AVLR[200][65];
	//unsigned long RecordID_AVLR[200],Reserved_AVRL[200];
	//long long RecordLenghtBeforeFooter_AVLR[200];

	//GeoKeyDirectory
	unsigned short KeyVersion[200],KeyRevision[200],MinorRevision[200],NumberofKeys[200];
	unsigned short TiffTagLocation[200],KeyCount[200], KeyOffset[200],KeyID[200]; // does not take 2d array

	//GeodoubleParams
	double geodouble[200];

	//GeoAsciiParams
	char geoascii[200];

	//SCANNER
	unsigned long SizeScanner[200],ReservedScanner[200],ScanPattern[200],NumberofMirrorFacets[200];
	unsigned char Instument[200][65], Serial[200][65]; // Allowing 100 more bytes in the scanner VLR (Scanner_additional)
	float Wavelength[200],PulseWidth[200],ScanFrequency[200],ScanAngleMin[200],ScanAngleMax[200],PulseFrequency[200],BeamDiameterAtExitAperture[200],BeamDivergence[200],MinimalRange[200],MaximalRange[200];
	char DescriptionScanner[200][65], Scanner_additional[200][101];

	//PULSE DESCRIPTOR
	//Composition record //200 size is the maximum number of AVLRs/VLRs
	unsigned long SizeComposition[200],Reserved_Composition[200],CompressionComposition[200],ScannerIndex[200];
	long OpticalCentreToAnchor[200];
	short NumberOfExtraWave[200],NumberOfSamplings[200],PulseDescrlocation;
	float SampleUnitsComposition[200];
	char DescriptionComposition[200][65],Discriptor_additional[200][100];// Allowing 100 more bytes in the table VLR (Discriptor_additional)

	//Sampling record // allowing 15 number of samplings
	unsigned long	SizeSampling[200][15],ReservedSampling[200][15],NumberOfSamples[200][15],CompressionSampling[200][15];
	unsigned char Type[200][15], Channel[200][15], Unused[200][15], BitsForDurationFromAnchor[200][15], BitsForNumberOfSegments[200][15],BitsForNumberOfSamples[200][15];
	float ScaleForDurationFromAnchor[200][15],OffsetForDurationFromAnchor[200][15],SampleUnits[200][15];
	char DescriptionSampling[200][15][65],Sampling_additional[200][100];
	unsigned short NumberOfSegments[200][15],BitsPerSample[200][15],LookupTableIndex[200][15];

	//TABLE //200 size is the maximum number of AVLRs/VLRs
	unsigned long SizeTableHR[200],ReservedTableHR[200],NumberOfTables[200],SizeTable[200][20],ReservedTable[200][20],NumberEntries[200][20],CompressionTable[200][20];// allowing 20 lookup tables
	unsigned char DataType[200][20], OptionsTable[200][20];// Allowing 100 more bytes in the table VLR (Table_additional)// allowing 20 lookup tables
	char DescriptionTable[200][20][65], DescriptionTableHR[200][65],TableHR_additional[200][100],Table_additional[200][100];// allowing 20 lookup tables
	unsigned short UnitOfMeasurement[200][20];// allowing 20 lookup tables

	  //WVS file header
	//char	FileSgnWvs[17];
	//unsigned long CompressionWvs;
	unsigned char ExtraWavesBytes[100];//allowing 100 extra bytes 
	//unsigned char ReservedWvs[45];//allowing 100 extra bytes 

	///Samplings in WVS // allowing 15 samplings, 15 segments and 300 samples
	//unsigned char Samples[15][15][300],DurationFromAnchor;
  	unsigned short NumberOfSegmentsWws[15];
	long NumberOfSamplesWvs[15] ;
	unsigned long Samples[15][1000];
	float DurationFromAnchorWvs[15],Offset_segment[15];
	
	int i,j,k,m; /* Loop index */
	double max_waveform_elevation;

	double peak_e, peak_n, peak_h; /* For point-cloud output */
	double sample_e, sample_n, sample_h; /* For wavecloud output */

	int waveform_count;
	int process_all; /* 1 = don't subset, 0 do subset */

	char input_file_param[2];
	int no_paramfile; /* For more than one pass within a swathe */

	double aircraft_altitude, aircraft_altitude_sum;
	double first_gps_time, last_gps_time;

	int bin_truncation_limit, peak_search_start;
	int waveform_bin_stats[5001]; /* Stats on how many bins in waveform */

	
	/* Mapping variables */
	int count_transmittedonly,count_multisegment,no_transmitted,low_saturated,count_disregarded,segment_too_short;

	/* Flightline descriptor variables */
	double westmost_aircraft_e, westmost_aircraft_n, eastmost_aircraft_e, eastmost_aircraft_n;

	//////////////////////////////////////////////////////////////////
	//unsigned long long int waveform_location;
	long long previous_gps_time;
	double start_pulse_easting, start_pulse_northing;
	double start_pulse_elevation, waveform_vector_easting, waveform_vector_northing, waveform_vector_elevation;
	unsigned short int waveform_length; 
	long start_pulse_length;
	char sample_depth; 

	double waveform_start_easting, waveform_start_northing, waveform_start_elevation;
	double waveform_end_easting, waveform_end_northing, waveform_end_elevation;

	int start_pulse_waveform[1001];
	int surface_return_waveform[9001], max_surface_return,max_start_pulse_return;
	//unsigned char waveform_sample_byte;
	int max_surface_bin,max_start_pulse_bin;
  
	//
	int max_surface_return_Sg[15],max_start_pulse_return_Sg[15],start_pulse_length_Sg[15];
	int max_surface_bin_Sg[15],max_start_pulse_bin_Sg[15],NumberOfSamplesWvs_Sg[15];
	float DurationFromAnchorWvs_Sg[15];
	int LastSample1;
	unsigned short **start_pulse_waveform_Sg;
	unsigned short **surface_return_waveform_Sg;

	surface_return_waveform_Sg=(unsigned short **)calloc(100,sizeof(unsigned short*));
	for (i=0;i<100;i++)
	{
		surface_return_waveform_Sg[i]=(unsigned short *) calloc(9000,sizeof(unsigned short));
	}

	start_pulse_waveform_Sg=(unsigned short **)calloc(100,sizeof(unsigned short*));
	for (i=0;i<100;i++)
	{
		start_pulse_waveform_Sg[i]=(unsigned short *) calloc(9000,sizeof(unsigned short));
	}

 
	count_disregarded=0;

  /* Main program */
	printf("***********************************************************************\n");
	printf("**               Module 1 of Canopy Height Program                   **\n");
	printf("***********************************************************************\n");
	printf("** Reads PulseWaves data format and converts to a .txt file          **\n");
	printf("** for further use with remaining modules of the program.            **\n");
	printf("**                                                                   **\n");
	printf("** Authors: Karolina Fieber, Ian Davenport                           **\n");
	printf("** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	printf("** University of Reading, UK                              March 2015 **\n");
	printf("***********************************************************************\n\n");



	printf("Enter the name of the input file (without extension): ");
	scanf("%s",&filename);

	sprintf(log_filename,"%s_info_Module1.txt",filename);
	log_fp=fopen(log_filename,"w");
	if (!log_fp)
	{
		printf("Can't open output info data file <<%s>>\n",log_filename);
		return(-1);
	}

	fprintf(log_fp,"***********************************************************************\n");
	fprintf(log_fp,"**               Module 1 of Canopy Height Program                   **\n");
	fprintf(log_fp,"***********************************************************************\n");
	fprintf(log_fp,"** Reads PulseWaves data format and converts to a .txt file          **\n");
	fprintf(log_fp,"** for further use with remaining modules of the program             **\n");
	fprintf(log_fp,"**                                                                   **\n");
	fprintf(log_fp,"** Authors: Karolina Fieber, Ian Davenport                           **\n");
	fprintf(log_fp,"** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	fprintf(log_fp,"** University of Reading, UK                              March 2015 **\n");
	fprintf(log_fp,"***********************************************************************\n\n");

	
	sprintf(pls_filename,"%s.pls",filename);
	sprintf(wvs_filename,"%s.wvs",filename);
	
	pls_fp=fopen(pls_filename,"rb");
	if (!pls_fp)
	{
		printf("Can't open input pls data file <<%s>>\n",pls_filename);
		fprintf(log_fp,"\nCan't open input pls data file <%s>\n",pls_filename);
		return(-1);
	}
	else
	{
		printf("\n\tReading from file <<%s>>\n",pls_filename);
		fprintf(log_fp,"\tReading from file <<%s>>\n",pls_filename);
	}

	wvs_fp=fopen(wvs_filename,"rb");
	if (!wvs_fp)
	{
		printf("Can't open input pls data file <<%s>>\n",wvs_filename);
		fprintf(log_fp,"Can't open input pls data file <<%s>>\n",wvs_filename);
		return(-1);
	}
	else
	{
		printf("\tReading from file <<%s>>\n",wvs_filename);
		fprintf(log_fp,"\tReading from file <<%s>>\n",wvs_filename);
	}


  //*********************************************************************************************************************************************************************************************//
  //*********************************************************************************************************************************************************************************************//
  //													// GENERATE EXTENT FILE //
  //*********************************************************************************************************************************************************************************************//
  //*********************************************************************************************************************************************************************************************//

	printf("\n***********************************************************************");
	printf("\n***********************  FIND EXTENTS OF THE DATA  ********************");
	printf("\n***********************************************************************\n\n");

	fprintf(log_fp,"\n***********************************************************************");
	fprintf(log_fp,"\n***********************  FIND EXTENTS OF THE DATA  ********************");
	fprintf(log_fp,"\n***********************************************************************\n\n");


	/* Scan through */

	min_easting = max_easting = min_northing = max_northing = first_gps_time = last_gps_time = min_elevation= max_elevation= -999 ;
	waveform_count=0;
	for (i=1;i<=5000;i++)
		waveform_bin_stats[i]=0;
	aircraft_altitude_sum=0;

	/////////////////////////////////////
	// READ HEADER OF THE DATA FILE
	////////////////////////////////////
	for (i=1;i<=16; i++)
		FileSgn[i]=readchar(pls_fp);
	GlobalParam=readulong(pls_fp);
	FileSourceID=readulong(pls_fp);
	ProjectID1=readulong(pls_fp);
	ProjectID2=readushort(pls_fp);
	ProjectID3=readushort(pls_fp);
	for (i=1;i<=8; i++)
		ProjectID4[i]=readuchar(pls_fp);
	for (i=1;i<=64; i++)
		SystemID[i]=readchar(pls_fp);
	for (i=1;i<=64; i++)
		Software[i]=readchar(pls_fp);
	CreationDayofYear=readushort(pls_fp);
	FileCreationYear=readushort(pls_fp);
	VersionMajor=readuchar(pls_fp);
	VersionMinor=readuchar(pls_fp);
	HeaderSize=readushort(pls_fp);
	OffsetToPulse=readlonglong(pls_fp);
	NumberOfPulses=readlonglong(pls_fp);
	PulseFormat=readulong(pls_fp);
	PulseAttributes=readulong(pls_fp);
	PulseSize=readulong(pls_fp);
	PulseCompression=readulong(pls_fp);
	ReservedPLS=readlonglong(pls_fp);
	NumberOfVLR=readulong(pls_fp);
	NumberOfAVLR=readlong(pls_fp);
	TScaleFactor=readdouble(pls_fp);
	TOffset=readdouble(pls_fp);
	MinT=readlonglong(pls_fp);
	MaxT=readlonglong(pls_fp);
	XScale=readdouble(pls_fp);
	YScale=readdouble(pls_fp);
	ZScale=readdouble(pls_fp);
	XOffset=readdouble(pls_fp);
	YOffset=readdouble(pls_fp);
	ZOffset=readdouble(pls_fp);
	MinX=readdouble(pls_fp);
	MaxX=readdouble(pls_fp);
	MinY=readdouble(pls_fp);
	MaxY=readdouble(pls_fp);
	MinZ=readdouble(pls_fp);
	MaxZ=readdouble(pls_fp);

	/////////////////////////////////////////
	// Read Variable Length Record (VLRs)
	/////////////////////////////////////////
	for (i=1; i<=(int)NumberOfVLR; i++)
	{
		//if(RecordID_VLR[i-1]==100001)
			//New_reserved[i]=readulong(pls_fp); // this field is not described in specification; before the first pulse descriptor record or after last scanner record
		for (j=1;j<=16; j++)
			UserID_VLR[i][j]=readchar(pls_fp);
		RecordID_VLR[i]=readulong(pls_fp);
		Reserved_VRL[i]=readulong(pls_fp);
		RecordLenghtAfterHeader_VLR[i]=readlonglong(pls_fp);
		for(j=1;j<=64; j++)
			Description_VLR[i][j]=readchar(pls_fp);

		//GeoKeyDirectory
		if (RecordID_VLR[i]==34735)
		{
			KeyVersion[i]=readushort(pls_fp);
			KeyRevision[i]=readushort(pls_fp);
			MinorRevision[i]=readushort(pls_fp);
			NumberofKeys[i]=readushort(pls_fp);
			for (j=1;j<=NumberofKeys[i];j++)
			{
				KeyID[i]=readushort(pls_fp);
				TiffTagLocation[i]=readushort(pls_fp);
				KeyCount[i]=readushort(pls_fp);
				KeyOffset[i]=readushort(pls_fp);
			}
		}
				
		//GeoDoubleParams
		if (RecordID_VLR[i]==34736)
		{
			for (j=1;j<=RecordLenghtAfterHeader_VLR[i]/8;j++)
			{
				geodouble[i]=readdouble(pls_fp);
			}
		}
					
		//GeoAsciiParams
		if (RecordID_VLR[i]==34737)
		{
			for (j=1;j<=RecordLenghtAfterHeader_VLR[i];j++)
			{
				geoascii[i]=readuchar(pls_fp);
			}
		}

					
		//SCANNER
		if (RecordID_VLR[i]>=100001&&RecordID_VLR[i]<100255) 
		{
			SizeScanner[i]=readulong(pls_fp);
			ReservedScanner[i]=readulong(pls_fp);
			for (j=1;j<=64; j++)
				Instument[i][j]=readuchar(pls_fp);
			for (j=1;j<=64; j++)
				Serial[i][j]=readuchar(pls_fp);
			Wavelength[i]=readfloat(pls_fp);
			PulseWidth[i]=readfloat(pls_fp);
			ScanPattern[i]=readulong(pls_fp);
			NumberofMirrorFacets[i]=readulong(pls_fp);
			ScanFrequency[i]=readfloat(pls_fp);
			ScanAngleMin[i]=readfloat(pls_fp);
			ScanAngleMax[i]=readfloat(pls_fp);
			PulseFrequency[i]=readfloat(pls_fp);
			BeamDiameterAtExitAperture[i]=readfloat(pls_fp);
			BeamDivergence[i]=readfloat(pls_fp);
			MinimalRange[i]=readfloat(pls_fp);
			MaximalRange[i]=readfloat(pls_fp);
			if (RecordLenghtAfterHeader_VLR[i]>SizeScanner[i])
			{
				for (j=1;j<=RecordLenghtAfterHeader_VLR[i]-248;j++) //total of 248 bytes in the default
					Scanner_additional[i][j]=readchar(pls_fp);
			}
			for (j=1;j<=64;j++)
				DescriptionScanner[i][j]=readchar(pls_fp);
		}

		//PULSE DESCRIPTOR
		if (RecordID_VLR[i]>=200001&&RecordID_VLR[i]<200255) 
		{
			//Composition record
			SizeComposition[i]=readulong(pls_fp);
			Reserved_Composition[i]=readulong(pls_fp);
			OpticalCentreToAnchor[i]=readlong(pls_fp);
			NumberOfExtraWave[i]=readushort(pls_fp);
			NumberOfSamplings[i]=readushort(pls_fp);
			SampleUnitsComposition[i]=readfloat(pls_fp);
			CompressionComposition[i]=readulong(pls_fp);
			ScannerIndex[i]=readulong(pls_fp);
			if ((int)SizeComposition[i]>92) // total of 92 bytes in composition record
			{
				for (j=1;j<=(int)SizeComposition[i]-92;j++)  
					Discriptor_additional[i][j]=readchar(pls_fp);
			}						
			//...
			for (j=1;j<=64;j++)	
				DescriptionComposition[i][j]=readchar(pls_fp); 

			//Sampling record
			for (j=1;j<=NumberOfSamplings[i];j++)
			{
				SizeSampling[i][j]=readulong(pls_fp);
				ReservedSampling[i][j]=readulong(pls_fp);
				Type[i][j]=readuchar(pls_fp); // 1 - outgoing, 2-returning
				Channel[i][j]=readuchar(pls_fp);
				Unused[i][j]=readuchar(pls_fp);
				BitsForDurationFromAnchor[i][j]=readuchar(pls_fp);
				ScaleForDurationFromAnchor[i][j]=readfloat(pls_fp);
				OffsetForDurationFromAnchor[i][j]=readfloat(pls_fp);
				BitsForNumberOfSegments[i][j]=readuchar(pls_fp);
				BitsForNumberOfSamples[i][j]=readuchar(pls_fp);
				NumberOfSegments[i][j]=readushort(pls_fp);
				NumberOfSamples[i][j]=readulong(pls_fp);
				BitsPerSample[i][j]=readushort(pls_fp);
				LookupTableIndex[i][j]=readushort(pls_fp);
				SampleUnits[i][j]=readfloat(pls_fp);
				CompressionSampling[i][j]=readulong(pls_fp);
				if ((int)SizeSampling[i][j]>104) //...
				{
					for (j=1;j<=(int)SizeSampling[i][j]-104;j++)  // total of 104 bytes in sampling records
						Sampling_additional[i][j]=readchar(pls_fp);
				}	
				for (k=1;k<=64;k++)	
					DescriptionSampling[i][j][k]=readchar(pls_fp); 
			}
			if (RecordID_VLR[i]==200001)
				PulseDescrlocation=i-1;
		}

		//TABLE
		if (RecordID_VLR[i]>=300001&&RecordID_VLR[i]<300255) 
		{
			SizeTableHR[i]=readulong(pls_fp); //header: 76 bytes in total
			ReservedTableHR[i]=readulong(pls_fp);
			NumberOfTables[i]=readulong(pls_fp);
			if ((int)SizeTableHR[i]>76)//...
			{
				for (j=1;j<=(int)SizeTableHR[i]-76;j++) // total of 76 bytes in default Table header
					TableHR_additional[i][j]=readchar(pls_fp);
			}	
			for (j=1;j<=64;j++)	
				DescriptionTableHR[i][j]=readchar(pls_fp);	

			for (j=1;j<=(int)NumberOfTables[i];j++)
			{
				SizeTable[i][j]=readulong(pls_fp); // total 84 bytes
				ReservedTable[i][j]=readulong(pls_fp);
				NumberEntries[i][j]= readulong(pls_fp);
				UnitOfMeasurement[i][j]=readushort(pls_fp);
				DataType[i][j]=readuchar(pls_fp);
				OptionsTable[i][j]=readuchar(pls_fp);
				CompressionTable[i][j]=readulong(pls_fp);
				if ((int)SizeTable[i][j]>84)//...
				{
					for (j=1;j<=(int)SizeTable[i][j]-84;j++) // total of 84 bytes in default Table header
						Table_additional[i][j]=readchar(pls_fp);
				}
				for(k=1;k<=64;k++)	
					DescriptionTable[i][j][k]=readchar(pls_fp);	
			}
		} // end of TABLE
	} //end of variable length records
	/////////////////////////////
	// Read Pulse Records
	/////////////////////////////
	for (i=1;i<=NumberOfPulses;i++)
	{
		GPSTime=readlonglong(pls_fp); /* T - GPS time tag of laser shot [GPS seconds of the week] */
		OffsetToWaves=readlonglong(pls_fp);
		AnchorX=readlong(pls_fp);
		AnchorY=readlong(pls_fp);
		AnchorZ=readlong(pls_fp);
		TargetX=readlong(pls_fp);
		TargetY=readlong(pls_fp);
		TargetZ=readlong(pls_fp);
		FirstSample=readshort(pls_fp);
		LastSample=readshort(pls_fp);
		PulseDescrIndex=readuchar(pls_fp);
		Reserved_pulse=readuchar(pls_fp);
		Intensity=readuchar(pls_fp);
		Classification=readuchar(pls_fp);
				  				  
		start_pulse_easting=AnchorX*XScale+XOffset; /* E0 */
		start_pulse_northing=AnchorY*YScale+YOffset ; /* N0 */
    	start_pulse_elevation=AnchorZ*ZScale+ZOffset;/* H0 */

		waveform_vector_easting=-(AnchorX-TargetX)*XScale/1000; /* dE */
		waveform_vector_northing=-(AnchorY-TargetY)*YScale/1000; /* dN */
		waveform_vector_elevation=-(AnchorZ-TargetZ)*ZScale/1000; /* dH */

		// E,N,Z coordinates of the begining and end of the waveform
		waveform_start_easting=start_pulse_easting+waveform_vector_easting*FirstSample;
    	waveform_start_northing=start_pulse_northing+waveform_vector_northing*FirstSample;
		waveform_start_elevation=start_pulse_elevation+waveform_vector_elevation*FirstSample;
				
		waveform_length=LastSample-FirstSample;
		waveform_end_easting=waveform_start_easting+waveform_vector_easting*waveform_length;
		waveform_end_northing=waveform_start_northing+waveform_vector_northing*waveform_length;
		waveform_end_elevation=waveform_start_elevation+waveform_vector_elevation*waveform_length;
				
		aircraft_altitude_sum+=start_pulse_elevation;

		if (waveform_start_northing>1e10)
			printf("Extremely large northing coordinate!\n");
		if ((waveform_start_easting<min_easting)||(min_easting<-990))
			min_easting=waveform_start_easting;
		if ((waveform_start_easting>max_easting)||(max_easting<-990))
			max_easting=waveform_start_easting;
		if ((waveform_start_northing<min_northing)||(min_northing<-990))
			min_northing=waveform_start_northing;
		if ((waveform_start_northing>max_northing)||(max_northing<-990))
			max_northing=waveform_start_northing;
		if ((waveform_end_easting<min_easting)||(min_easting<-990))
			min_easting=waveform_end_easting;
		if ((waveform_end_easting>max_easting)||(max_easting<-990))
			max_easting=waveform_end_easting;
		if ((waveform_end_northing<min_northing)||(min_northing<-990))
			min_northing=waveform_end_northing;
		if ((waveform_end_northing>max_northing)||(max_northing<-990))
			max_northing=waveform_end_northing;
		if ((GPSTime<first_gps_time)||(first_gps_time<-990))
			first_gps_time=GPSTime;
		if ((GPSTime>last_gps_time)||(last_gps_time<-990))
			last_gps_time=GPSTime;
		if (max_elevation<waveform_start_elevation||(max_elevation<-990))
			max_elevation=waveform_start_elevation;
		if (min_elevation>waveform_end_elevation||(min_elevation<-990))
			min_elevation=waveform_end_elevation;
	
		if (waveform_length<=5000)
			waveform_bin_stats[(int)waveform_length]++;
		else
			printf("Waveform > 5000 bins!\n");
		waveform_count++;		
	} // end of pulse records


	aircraft_altitude=aircraft_altitude_sum/waveform_count;
	printf("* Number of waveforms found: \t%d *\n",waveform_count);
	//printf("* Extent: \t\t\t%.3lf - %.3lf E, \n\t\t\t\t%.3lf - %.3lf N * \n* Mean aircraft elevation: \t%.2lfm * \n* GPS time: \t\t\t%lf-%lf *\n",min_northing, max_northing, min_easting, max_easting,aircraft_altitude,first_gps_time, last_gps_time);
	printf("* Extent Easting: \t\t%.3lf - %.3lf E *\n* Extent Northing:\t\t%.3lf - %.3lf N *",min_easting, max_easting,min_northing, max_northing);
	printf("\n* Extent Elevation: \t\t%.3lf - %.3lf *",min_elevation, max_elevation);
	printf("\n* Mean anchor elevation: \t%.3lfm * \n* GPS time: \t\t\t%.2lf-%.2lf *\n",aircraft_altitude,first_gps_time, last_gps_time);
	
	fprintf(log_fp,"DATA INFORMATION:");
	fprintf(log_fp,"\n* Number of waveforms found: \t%d *\n",waveform_count);
	fprintf(log_fp,"* Extent Easting: \t\t%.3lf - %.3lf E *\n* Extent Northing:\t\t%.3lf - %.3lf N *",min_easting, max_easting,min_northing, max_northing);
	fprintf(log_fp,"\n* Extent Elevation: \t\t%.3lf - %.3lf",min_elevation, max_elevation);
	fprintf(log_fp,"\n* Mean anchor elevation: \t%.3lfm * \n* GPS time: \t\t\t%.2lf-%.2lf *\n",aircraft_altitude,first_gps_time, last_gps_time);
	

   //*******************************************************//
				// END OF FIND EXTENTS OF THE DATA //
  //*******************************************************//

 
  //*********************************************************************************************************************************************************************************************//
  //*********************************************************************************************************************************************************************************************//
																// EXTRACT SELECTED AREA //
  //*********************************************************************************************************************************************************************************************//
  //*********************************************************************************************************************************************************************************************//
  
	/* Read area and output from parameter file */
    /* Two ways of specifying area - either min/max northing/easting or centre and radius */

	printf("\n***********************************************************************");
	printf("\n************************  EXTRACT SELECTED AREA  **********************");
	printf("\n***********************************************************************\n\n");

	fprintf(log_fp,"\n***********************************************************************");
	fprintf(log_fp,"\n************************  EXTRACT SELECTED AREA  **********************");
	fprintf(log_fp,"\n***********************************************************************\n\n");


	bin_truncation_limit=1000;  

	printf("Do you want to use %s_parameter_area input file (Y/N): ",filename); 
	fprintf(log_fp,"Do you want to use %s_parameter_area input file (Y/N): ",filename);
	scanf("%s",&input_file_param);
	while ((input_file_param[0]!='Y'&&input_file_param[0]!='y')&&(input_file_param[0]!='N'&&input_file_param[0]!='n'))
	{
		printf("\nWRONG INPUT: Y/N? ");
		scanf("%s",&input_file_param);
	}
	fprintf(log_fp,"%s\n",&input_file_param);

	if ((input_file_param[0]=='Y')||(input_file_param[0]=='y'))
	{
		sprintf(param_filename,"%s_parameter_area.txt",filename);
	}
	 
	no_paramfile=0;
	if ((input_file_param[0]=='Y')||(input_file_param[0]=='y'))
	{
		param_fp=fopen(param_filename,"r");
		if (!param_fp)
		{	  
			printf("\tCan't open parameter file <<%s>>\n",param_filename);
			fprintf(log_fp,"\tCan't open parameter file <<%s>>\n",param_filename);
			if (!param_fp)
			{
				no_paramfile=1;
			}
		}
		else
		{
			no_paramfile=0;
		}
	}
	else
	{
		no_paramfile=1;
	}

	if (no_paramfile==0) //if parameter file exists
	{
		fscanf(param_fp,"%d",&area_type); /* 0 = all, 1 = circle, 2 = rectangle, 3 = cuboid */
	  
		fprintf(log_fp, "Reading from <<%s>>\n",param_filename); 
		if (area_type==1) /* circle */
		{
			fscanf(param_fp,"%lf",&easting_centre);
			fscanf(param_fp,"%lf",&northing_centre);
			fscanf(param_fp,"%lf",&radius);
			printf("\nExtracting within %4.1lfm of %9.1lfN, %10.1lfE\n",radius,northing_centre,easting_centre);
		}
		else if (area_type==2) /* rectangle */
		{

			fscanf(param_fp,"%lf",&west_edge);
			fscanf(param_fp,"%lf",&east_edge);
			fscanf(param_fp,"%lf",&south_edge);
			fscanf(param_fp,"%lf",&north_edge);
			easting_centre=(west_edge+east_edge)/2;
			northing_centre=(south_edge+north_edge)/2;
		}
		else if (area_type==3) /* cuboid */
		{

			fscanf(param_fp,"%lf",&west_edge);
			fscanf(param_fp,"%lf",&east_edge);
			fscanf(param_fp,"%lf",&south_edge);
			fscanf(param_fp,"%lf",&north_edge);
			fscanf(param_fp,"%lf",&low_edge);
			fscanf(param_fp,"%lf",&high_edge);
			fscanf(param_fp,"%d",&extract_threshold); // minimum amplitude that a single bin has to have wihtin the cuboid in order for the entire waveform to be extracted
		}

		if (area_type==0)
			process_all=1;
		else
			process_all=0;
	}

	if (no_paramfile==1)
	{
		printf("\nSpecify area of interest: \t0 = all, \n\t\t\t\t1 = circle, \n\t\t\t\t2 = rectangle, \n\t\t\t\t3 = cuboid:  "); 
		scanf("%d",&area_type);
		fprintf(log_fp, "\nSpecify area of interest (0=all, 1=circle, 2=rectangle, 3=cuboid):  %d\n",area_type); 
		
		if (area_type==1) /* circle */
		{
			printf("\nSpecify Easting coordinate of the site's centre [m]:  "); 
			scanf("%lf",&easting_centre);
			printf("Specify Northing coordinate of the site's centre [m]:  "); 
			scanf("%lf",&northing_centre);
			printf("Specify radius of the site [m]:  "); 
			scanf("%lf",&radius);
			fprintf(log_fp,"Extracting within %4.1lfm of %10.1lfE and %9.1lfN\n",radius,easting_centre,northing_centre);
		}
		else if (area_type==2) /* rectangle */
		{
			printf("\nSpecify West coordinate of the site [m]:  "); 
			scanf("%lf",&west_edge);
			printf("Specify East coordinate of the site [m]:  "); 
			scanf("%lf",&east_edge);
			printf("Specify South coordinate of the site [m]:  "); 
			scanf("%lf",&south_edge);
			printf("Specify North coordinate of the site [m]:  "); 
			scanf("%lf",&north_edge);

			easting_centre=(west_edge+east_edge)/2;
			northing_centre=(south_edge+north_edge)/2;
			fprintf(log_fp,"Extracting within %9.1lfE and %10.1lfE and %9.1lfN and %10.1lfN \n",west_edge,east_edge,south_edge,north_edge);
		}
		else if (area_type==3) /* cuboid */
		{
			printf("\nSpecify West coordinate of the site [m]:  "); 
			scanf("%lf",&west_edge);
			printf("Specify East coordinate of the site [m]:  "); 
			scanf("%lf",&east_edge);
			printf("Specify South coordinate of the site [m]:  "); 
			scanf("%lf",&south_edge);
			printf("Specify North coordinate of the site [m]:  "); 
			scanf("%lf",&north_edge);
			printf("Specify Low boundary of the site [m]:  "); 
			scanf("%lf",&low_edge);
			printf("Specify High boundary of the site [m]:  "); 
			scanf("%lf",&high_edge);
			printf("Specify amplitude threshold - waveforms with at least one bin above that value will be extracted:  "); 
			scanf("%d",&extract_threshold);
			fprintf(log_fp,"Extracting within %9.1lfE and %10.1lfE and %9.1lfN and %10.1lfN \n",west_edge,east_edge,south_edge,north_edge);
			fprintf(log_fp,"Elevation boundary: %lf and %lf\n", low_edge, high_edge);
		}
		if (area_type==0)
			process_all=1;
		else
			process_all=0;
	}

	if (area_type!=1 && area_type!=2 && area_type!=3 && area_type!=0)
	{
		printf("WRONG input of the area type. Processing all.\n"); 
		fprintf(log_fp,"\nWRONG input of the area type. Processing all."); 
		process_all=1;
		//area_type=0;
	}
	

	/* First check if potential, then search */
	if (
		((easting_centre>min_easting)&&(easting_centre<max_easting)&&(northing_centre>min_northing)&&(northing_centre<max_northing))
		||
		/* SW */ ((area_type>1)&&(west_edge>min_easting)&&(west_edge<max_easting)&&(south_edge>min_northing)&&(south_edge<max_northing))
		||
		/* SE */ ((area_type>1)&&(east_edge>min_easting)&&(east_edge<max_easting)&&(south_edge>min_northing)&&(south_edge<max_northing))
		||
		/* NW */ ((area_type>1)&&(west_edge>min_easting)&&(west_edge<max_easting)&&(north_edge>min_northing)&&(north_edge<max_northing))
		||
		/* NE */ ((area_type>1)&&(east_edge>min_easting)&&(east_edge<max_easting)&&(north_edge>min_northing)&&(north_edge<max_northing))

		||

		/* SW */ ((area_type>1)&&(min_easting>west_edge)&&(min_easting<east_edge)&&(min_northing>south_edge)&&(min_northing<north_edge))
		||
		/* SE */ ((area_type>1)&&(max_easting>west_edge)&&(max_easting<east_edge)&&(min_northing>south_edge)&&(min_northing<north_edge))
		||
		/* NW */ ((area_type>1)&&(min_easting>west_edge)&&(min_easting<east_edge)&&(max_northing>south_edge)&&(max_northing<north_edge))
		||
		/* NE */ ((area_type>1)&&(max_easting>west_edge)&&(max_easting<east_edge)&&(max_northing>south_edge)&&(max_northing<north_edge))
		|| (process_all==1)
		)
	{
		pls_fp=fopen(pls_filename,"rb");
	    if (!pls_fp)
	    {
		    printf("Can't open input pls data file <<%s>>\n",pls_filename);
			fprintf(log_fp,"\nCan't open input pls data file <<%s>>\n",pls_filename);
		    return(-1);
	    }


    	wvs_fp=fopen(wvs_filename,"rb");
	    if (!wvs_fp)
		{
			printf("Can't open input wvs data file <<%s>>\n",wvs_filename);
			fprintf(log_fp,"\nCan't open input wvs data file <<%s>>\n",wvs_filename);
		    return(-1);
		}


		sprintf(out_waveform_filename,"%s_fulwvs.txt",filename);
		out_waveform_fp=fopen(out_waveform_filename,"w");	
		if (!out_waveform_fp)
		{
			printf("Can't open output waveform data file <<%s>>\n",out_waveform_filename);
			fprintf(log_fp,"\nCan't open output waveform data file <<%s>>\n",out_waveform_filename);
			return(-1);
		}
			  
    	/* printf("Now reading and outputting waveforms...\n"); */

		waveform_count=0;
		for (i=1;i<=5000;i++)
			waveform_bin_stats[i]=0;
		/* Reset flight extremes */
		westmost_aircraft_e=-999;
		westmost_aircraft_n=-999;
		eastmost_aircraft_e=-999;
		eastmost_aircraft_n=-999;
			  
		/////////////////////////////////////
		// READ HEADER OF THE DATA FILE
		////////////////////////////////////
		for (i=1;i<=16; i++)
			FileSgn[i]=readchar(pls_fp);
		GlobalParam=readulong(pls_fp);
		FileSourceID=readulong(pls_fp);
		ProjectID1=readulong(pls_fp);
		ProjectID2=readushort(pls_fp);
		ProjectID3=readushort(pls_fp);
		for (i=1;i<=8; i++)
			ProjectID4[i]=readuchar(pls_fp);
		for (i=1;i<=64; i++)
			SystemID[i]=readchar(pls_fp);
		for (i=1;i<=64; i++)
			Software[i]=readchar(pls_fp);
		CreationDayofYear=readushort(pls_fp);
		FileCreationYear=readushort(pls_fp);
		VersionMajor=readuchar(pls_fp);
		VersionMinor=readuchar(pls_fp);
		HeaderSize=readushort(pls_fp);
		OffsetToPulse=readlonglong(pls_fp);
		NumberOfPulses=readlonglong(pls_fp);
		PulseFormat=readulong(pls_fp);
		PulseAttributes=readulong(pls_fp);
		PulseSize=readulong(pls_fp);
		PulseCompression=readulong(pls_fp);
		ReservedPLS=readlonglong(pls_fp);
		NumberOfVLR=readulong(pls_fp);
		NumberOfAVLR=readlong(pls_fp);
		TScaleFactor=readdouble(pls_fp);
		TOffset=readdouble(pls_fp);
		MinT=readlonglong(pls_fp);
		MaxT=readlonglong(pls_fp);
		XScale=readdouble(pls_fp);
		YScale=readdouble(pls_fp);
		ZScale=readdouble(pls_fp);
		XOffset=readdouble(pls_fp);
		YOffset=readdouble(pls_fp);
		ZOffset=readdouble(pls_fp);
		MinX=readdouble(pls_fp);
		MaxX=readdouble(pls_fp);
		MinY=readdouble(pls_fp);
		MaxY=readdouble(pls_fp);
		MinZ=readdouble(pls_fp);
		MaxZ=readdouble(pls_fp);

		/////////////////////////////////////////
		// Read Variable Length Record (VLRs)
		/////////////////////////////////////////
		for (i=1; i<=(int)NumberOfVLR; i++)
		{
			//if(RecordID_VLR[i-1]==100001)
				//New_reserved[i]=readulong(pls_fp); // this field is not described in specification; before the first pulse descriptor record or after last scanner record
			for (j=1;j<=16; j++)
				UserID_VLR[i][j]=readchar(pls_fp);
			RecordID_VLR[i]=readulong(pls_fp);
			Reserved_VRL[i]=readulong(pls_fp);
			RecordLenghtAfterHeader_VLR[i]=readlonglong(pls_fp);
			for(j=1;j<=64; j++)
				Description_VLR[i][j]=readchar(pls_fp);

			//GeoKeyDirectory
			if (RecordID_VLR[i]==34735)
			{
				KeyVersion[i]=readushort(pls_fp);
				KeyRevision[i]=readushort(pls_fp);
				MinorRevision[i]=readushort(pls_fp);
				NumberofKeys[i]=readushort(pls_fp);
				for (j=1;j<=NumberofKeys[i];j++)
				{
					KeyID[i]=readushort(pls_fp);
					TiffTagLocation[i]=readushort(pls_fp);
					KeyCount[i]=readushort(pls_fp);
					KeyOffset[i]=readushort(pls_fp);
				}
			}

			//GeoDoubleParams
			if (RecordID_VLR[i]==34736)
			{
				for (j=1;j<=RecordLenghtAfterHeader_VLR[i]/8;j++)
				{
					geodouble[i]=readdouble(pls_fp);
				}
			}
					
			//GeoAsciiParams
			if (RecordID_VLR[i]==34737)
			{
				for (j=1;j<=RecordLenghtAfterHeader_VLR[i];j++)
				{
					geoascii[i]=readuchar(pls_fp);
				}
			}
					
			//SCANNER
			if (RecordID_VLR[i]>=100001&&RecordID_VLR[i]<100255) 
			{
				SizeScanner[i]=readulong(pls_fp);
				ReservedScanner[i]=readulong(pls_fp);
				for (j=1;j<=64; j++)
					Instument[i][j]=readuchar(pls_fp);
				for (j=1;j<=64; j++)
					Serial[i][j]=readuchar(pls_fp);
				Wavelength[i]=readfloat(pls_fp);
				PulseWidth[i]=readfloat(pls_fp);
				ScanPattern[i]=readulong(pls_fp);
				NumberofMirrorFacets[i]=readulong(pls_fp);
				ScanFrequency[i]=readfloat(pls_fp);
				ScanAngleMin[i]=readfloat(pls_fp);
				ScanAngleMax[i]=readfloat(pls_fp);
				PulseFrequency[i]=readfloat(pls_fp);
				BeamDiameterAtExitAperture[i]=readfloat(pls_fp);
				BeamDivergence[i]=readfloat(pls_fp);
				MinimalRange[i]=readfloat(pls_fp);
				MaximalRange[i]=readfloat(pls_fp);
				if (RecordLenghtAfterHeader_VLR[i]>SizeScanner[i])
				{
					for (j=1;j<=RecordLenghtAfterHeader_VLR[i]-248;j++) //total of 248 bytes in the default
						Scanner_additional[i][j]=readchar(pls_fp);
				}
				for (j=1;j<=64;j++)
					DescriptionScanner[i][j]=readchar(pls_fp);
			}

			//PULSE DESCRIPTOR
			if (RecordID_VLR[i]>=200001&&RecordID_VLR[i]<200255) 
			{
				//Composition record
				SizeComposition[i]=readulong(pls_fp);
				Reserved_Composition[i]=readulong(pls_fp);
				OpticalCentreToAnchor[i]=readlong(pls_fp);
				NumberOfExtraWave[i]=readushort(pls_fp);
				NumberOfSamplings[i]=readushort(pls_fp);
				SampleUnitsComposition[i]=readfloat(pls_fp);
				CompressionComposition[i]=readulong(pls_fp);
				ScannerIndex[i]=readulong(pls_fp);
				if (SizeComposition[i]>92) // total of 92 bytes in composition record
				{
					for (j=1;j<=(int)SizeComposition[i]-92;j++)  
						Discriptor_additional[i][j]=readchar(pls_fp);
				}						
				//...
				for (j=1;j<=64;j++)	
					DescriptionComposition[i][j]=readchar(pls_fp); 

				//Sampling record
				for (j=1;j<=NumberOfSamplings[i];j++)
				{
					SizeSampling[i][j]=readulong(pls_fp);
					ReservedSampling[i][j]=readulong(pls_fp);
					Type[i][j]=readuchar(pls_fp); // 1 - outgoing, 2-returning
					Channel[i][j]=readuchar(pls_fp);
					Unused[i][j]=readuchar(pls_fp);
					BitsForDurationFromAnchor[i][j]=readuchar(pls_fp);
					ScaleForDurationFromAnchor[i][j]=readfloat(pls_fp);
					OffsetForDurationFromAnchor[i][j]=readfloat(pls_fp);
					BitsForNumberOfSegments[i][j]=readuchar(pls_fp);
					BitsForNumberOfSamples[i][j]=readuchar(pls_fp);
					NumberOfSegments[i][j]=readushort(pls_fp);
					NumberOfSamples[i][j]=readulong(pls_fp);
					BitsPerSample[i][j]=readushort(pls_fp);
					LookupTableIndex[i][j]=readushort(pls_fp);
					SampleUnits[i][j]=readfloat(pls_fp);
					CompressionSampling[i][j]=readulong(pls_fp);
					if (SizeSampling[i][j]>104) //...
					{
						for (j=1;j<=(int)SizeSampling[i][j]-104;j++)  // total of 104 bytes in sampling records
							Sampling_additional[i][j]=readchar(pls_fp);
					}	
					for (k=1;k<=64;k++)	
						DescriptionSampling[i][j][k]=readchar(pls_fp); 
				}
				if (RecordID_VLR[i]==200001) // check how many VLRs are before the first pulse descriptor
					PulseDescrlocation=i-1;
			}

			//TABLE
			if (RecordID_VLR[i]>=300001&&RecordID_VLR[i]<300255) 
			{
				SizeTableHR[i]=readulong(pls_fp); //header: 76 bytes in total
				ReservedTableHR[i]=readulong(pls_fp);
				NumberOfTables[i]=readulong(pls_fp);
				if (SizeTableHR[i]>76)//...
				{
					for (j=1;j<=(int)SizeTableHR[i]-76;j++) // total of 76 bytes in default Table header
						TableHR_additional[i][j]=readchar(pls_fp);
				}	
				for (j=1;j<=64;j++)	
					DescriptionTableHR[i][j]=readchar(pls_fp);	

				for (j=1;j<=(int)NumberOfTables[i];j++)
				{
					SizeTable[i][j]=readulong(pls_fp); // total 84 bytes
					ReservedTable[i][j]=readulong(pls_fp);
					NumberEntries[i][j]= readulong(pls_fp);
					UnitOfMeasurement[i][j]=readushort(pls_fp);
					DataType[i][j]=readuchar(pls_fp);
					OptionsTable[i][j]=readuchar(pls_fp);
					CompressionTable[i][j]=readulong(pls_fp);
					if (SizeTable[i][j]>84)//...
					{
						for (j=1;j<=(int)SizeTable[i][j]-84;j++) // total of 84 bytes in default Table header
							Table_additional[i][j]=readchar(pls_fp);
					}
					for(k=1;k<=64;k++)	
						DescriptionTable[i][j][k]=readchar(pls_fp);	
				}
			} //end of table
		} // end of variable length record
		/////////////////////////////
		// Read Pulse Records
		/////////////////////////////
		count_transmittedonly=0; //waveforms with only transmitted pulse recorded (no returning waveform)
		count_multisegment=0; // waveforms with more than one returning segment
		for (i=1;i<=NumberOfPulses;i++)
		{
			GPSTime=readlonglong(pls_fp); /* T - GPS time tag of laser shot [GPS seconds of the week] */
			OffsetToWaves=readlonglong(pls_fp);
			AnchorX=readlong(pls_fp);
			AnchorY=readlong(pls_fp);
			AnchorZ=readlong(pls_fp);
			TargetX=readlong(pls_fp);
			TargetY=readlong(pls_fp);
			TargetZ=readlong(pls_fp);
			FirstSample=readshort(pls_fp);
			LastSample=readshort(pls_fp);
			LastSample1=LastSample;
			PulseDescrIndex=readuchar(pls_fp);

			Reserved_pulse=readuchar(pls_fp);
			Intensity=readuchar(pls_fp);
			Classification=readuchar(pls_fp);
				  				  
			start_pulse_easting=AnchorX*XScale+XOffset; /* E0 */
			start_pulse_northing=AnchorY*YScale+YOffset ; /* N0 */
    		start_pulse_elevation=(AnchorZ*ZScale+ZOffset);/* H0 */

			waveform_vector_easting=-((AnchorX-TargetX)*XScale/1000); /* dE */
		    waveform_vector_northing=-((AnchorY-TargetY)*YScale/1000); /* dN */
			waveform_vector_elevation=-((AnchorZ-TargetZ)*ZScale/1000); /* dH */


			waveform_length=LastSample-FirstSample;
				 
			////////////////////////////////////////////
		    	//READ THE WVS FILE
			////////////////////////////////////////////
			/* Read in whole waveform, calc peak here */
	    	_fseeki64 (wvs_fp , (OffsetToWaves) , SEEK_SET );  /* Untested code, for files >4Gb */ // moves the pointer to the specific location in the file


			if (waveform_length>1000)
				peak_search_start=waveform_length/2;
			else
				peak_search_start=1;

			// Waves trains
			if (NumberOfExtraWave[PulseDescrIndex+PulseDescrlocation]!=0) // empty at the moment
			{
				for (j=1;j<=NumberOfExtraWave[PulseDescrIndex+PulseDescrlocation];j++) /* there was -1 after PulseDescrlocation not sure why*/
					ExtraWavesBytes[j]=readuchar(wvs_fp);
			}
					
			// Number of Samplings Loop
			for (j=1;j<=NumberOfSamplings[PulseDescrIndex+PulseDescrlocation]; j++) 
			{
				if (BitsForNumberOfSegments[PulseDescrIndex+PulseDescrlocation][j]!=0)
				{	
					if (BitsForNumberOfSegments[PulseDescrIndex+PulseDescrlocation][j]==8)
						NumberOfSegmentsWws[j]=readuchar(wvs_fp);
					if (BitsForNumberOfSegments[PulseDescrIndex+PulseDescrlocation][j]==16)
						NumberOfSegmentsWws[j]=readushort(wvs_fp);
				}
				else /* if BitsForNumberOfSegments[PulseDescrIndex+PulseDescrlocation][j]==0, segmenting is fixxed and specified by nymber of segments in sampling records*/
				{
					NumberOfSegmentsWws[j]=NumberOfSegments[PulseDescrIndex+PulseDescrlocation][j];
				}

				for (m=1;m<=NumberOfSegmentsWws[j]; m++) 
				{
					max_surface_bin=0;
    	    		max_surface_return=0;
					max_start_pulse_return=0;
					max_start_pulse_bin=0;
					//offset from the origin
					if (BitsForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]!=0)
					{
						if (BitsForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]==8) 
						{	
							DurationFromAnchorWvs[j]=readuchar(wvs_fp);
							DurationFromAnchorWvs[j]=DurationFromAnchorWvs[j]*ScaleForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]+OffsetForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j];
						}
						if (BitsForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]==16)
						{
							DurationFromAnchorWvs[j]=readushort(wvs_fp);
							DurationFromAnchorWvs[j]=DurationFromAnchorWvs[j]*ScaleForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]+OffsetForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j];
						}
						if (BitsForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]==32)
						{
							DurationFromAnchorWvs[j]=(float)readulong(wvs_fp);
							DurationFromAnchorWvs[j]=DurationFromAnchorWvs[j]*ScaleForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]+OffsetForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j];
						}
					}
					else if (BitsForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]==0) // recording from the place of origin
					{
						DurationFromAnchorWvs[j]=0*ScaleForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j]+OffsetForDurationFromAnchor[PulseDescrIndex+PulseDescrlocation][j];
					}
					DurationFromAnchorWvs_Sg[m]=DurationFromAnchorWvs[j]; //

					// VARIABLE SAMPLING CASE
					if (BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]!=0)//(BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]/BitsPerSample[PulseDescrIndex+PulseDescrlocation][j]!=0)
					{	
						if (BitsPerSample[PulseDescrIndex+PulseDescrlocation][j]==8) // one byte samples
						{
							if (BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]==8)
								NumberOfSamplesWvs[j]=readuchar(wvs_fp); 
							if (BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]==16)
								NumberOfSamplesWvs[j]=readshort(wvs_fp); 
							sample_depth=0;

							NumberOfSamplesWvs_Sg[m]=NumberOfSamplesWvs[j]; //

							for (k=1;k<=(int)NumberOfSamplesWvs[j];k++)
							{
								Samples[j][k]=readuchar(wvs_fp);
								if (Type[PulseDescrIndex+PulseDescrlocation][j]==1) //if the outgoing waveform
								{
									start_pulse_waveform[k]=(unsigned int)(Samples[j][k]);
									start_pulse_waveform_Sg[m][k]=(unsigned int)(Samples[j][k]);//
									start_pulse_length=NumberOfSamplesWvs[j];
									start_pulse_length_Sg[m]=start_pulse_length;
									if ((start_pulse_waveform[k]>max_start_pulse_return))
            						{
										max_start_pulse_return=start_pulse_waveform[k];
										max_start_pulse_return_Sg[m]=start_pulse_waveform[k]; //
		            					max_start_pulse_bin_Sg[m]=k; //
									}
								}
								else                                                // if returning waveform
								{
									surface_return_waveform[k]=(int)Samples[j][k];
									surface_return_waveform_Sg[m][k]=(int)Samples[j][k]; //
									if ((surface_return_waveform[k]>max_surface_return)&&(k>=peak_search_start))
            						{
										max_surface_return=surface_return_waveform[k];
										max_surface_return_Sg[m]=surface_return_waveform[k]; //
		            					max_surface_bin_Sg[m]=k;		//
									}
								}
							}
						}
						if (BitsPerSample[PulseDescrIndex+PulseDescrlocation][j]==16) // two byte samples - only returning waveform
						{
							if (BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]==8)
								NumberOfSamplesWvs[j]=readuchar(wvs_fp); 
							if (BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]==16)
								NumberOfSamplesWvs[j]=readshort(wvs_fp); 
							sample_depth=1;
									
							NumberOfSamplesWvs_Sg[m]=NumberOfSamplesWvs[j]; //

							for (k=1;k<=(int)NumberOfSamplesWvs[j];k++)
							{
								Samples[j][k]=readushort(wvs_fp);
								if (Type[PulseDescrIndex+PulseDescrlocation][j]==1) //if the outgoing waveform
								{
									start_pulse_length=NumberOfSamplesWvs[j];
									start_pulse_length_Sg[m]=start_pulse_length;
									start_pulse_waveform[k]=(unsigned int)(Samples[j][k]);
									start_pulse_waveform_Sg[m][k]=(unsigned int)(Samples[j][k]);//
									if ((start_pulse_waveform[k]>max_start_pulse_return))
            						{
										max_start_pulse_return=start_pulse_waveform[k];
										max_start_pulse_return_Sg[m]=start_pulse_waveform[k]; //
		            					max_start_pulse_bin_Sg[m]=k; //
									}
								}
								else                                                // if returning waveform
								{
									surface_return_waveform[k]=(int)Samples[j][k];
									surface_return_waveform_Sg[m][k]=(int)Samples[j][k]; //
									if ((surface_return_waveform[k]>max_surface_return)&&(k>=peak_search_start))
            						{
										max_surface_return=surface_return_waveform[k];
										max_surface_return_Sg[m]=surface_return_waveform[k];
		            					max_surface_bin_Sg[m]=k;
									}
								}
							}
						}
					}
					//FIXED SAMPLING CASE
					else if (BitsForNumberOfSamples[PulseDescrIndex+PulseDescrlocation][j]==0)
					{
						if (BitsPerSample[PulseDescrIndex+PulseDescrlocation][j]==8) // one byte samples
						{
							NumberOfSamplesWvs[j]=NumberOfSamples[PulseDescrIndex+PulseDescrlocation][j];
							sample_depth=0;

							NumberOfSamplesWvs_Sg[m]=NumberOfSamplesWvs[j]; //

							for (k=1;k<=(int)NumberOfSamplesWvs[j];k++)
							{
								Samples[j][k]=readuchar(wvs_fp);
								if (Type[PulseDescrIndex+PulseDescrlocation][j]==1) //if the outgoing waveform
								{
									start_pulse_length=NumberOfSamplesWvs[j];
									start_pulse_length_Sg[m]=start_pulse_length;
									start_pulse_waveform[k]=(unsigned int)(Samples[j][k]);
									start_pulse_waveform_Sg[m][k]=(unsigned int)(Samples[j][k]);//
									if ((start_pulse_waveform[k]>max_start_pulse_return))
            						{
										max_start_pulse_return=start_pulse_waveform[k];
										max_start_pulse_return_Sg[m]=start_pulse_waveform[k];//
		            					max_start_pulse_bin_Sg[m]=k;//
									}
								}
								else                                                // if returning waveform
								{
									surface_return_waveform[k]=(int)Samples[j][k];
									surface_return_waveform_Sg[m][k]=(int)Samples[j][k]; //
									if ((surface_return_waveform[k]>max_surface_return)&&(k>=peak_search_start))
            						{
										max_surface_return=surface_return_waveform[k];
										max_surface_return_Sg[m]=surface_return_waveform[k]; //
		            					max_surface_bin_Sg[m]=k; //
									}
								}
							}
						}
						if (BitsPerSample[PulseDescrIndex+PulseDescrlocation][j]==16) // two byte samples - only returning waveform
						{
							NumberOfSamplesWvs[j]=NumberOfSamples[PulseDescrIndex+PulseDescrlocation][j];
							sample_depth=1;
							NumberOfSamplesWvs_Sg[m]=NumberOfSamplesWvs[j]; //

							for (k=1;k<=(int)NumberOfSamplesWvs[j];k++)
							{
								Samples[j][k]=readushort(wvs_fp);
										
								if (Type[PulseDescrIndex+PulseDescrlocation][j]==1) //if the outgoing waveform
								{
									start_pulse_length=NumberOfSamplesWvs[j];
									start_pulse_length_Sg[m]=start_pulse_length;
									start_pulse_waveform[k]=(unsigned int)(Samples[j][k]);
									start_pulse_waveform_Sg[m][k]=(unsigned int)(Samples[j][k]);//
									if ((start_pulse_waveform[k]>max_start_pulse_return))
            						{
										max_start_pulse_return=start_pulse_waveform[k];
										max_start_pulse_return_Sg[m]=start_pulse_waveform[k];//
		            					max_start_pulse_bin_Sg[m]=k;//
									}
								}
								else                                                // if returning waveform
								{
									surface_return_waveform[k]=(int)Samples[j][k];
									surface_return_waveform_Sg[m][k]=(int)Samples[j][k]; //
									if ((surface_return_waveform[k]>max_surface_return)&&(k>=peak_search_start))
            						{
										max_surface_return=surface_return_waveform[k];
										max_surface_return_Sg[m]=surface_return_waveform[k];//
		            					max_surface_bin_Sg[m]=k;//
									}
								}
							}
						}
					}
				}
				//for multiple segments DurationFromAnchor is sometimes missing the offset value
				if (NumberOfSegmentsWws[j]>1)
				{						
					if (DurationFromAnchorWvs_Sg[1]<FirstSample)
					{
						Offset_segment[1]=FirstSample-DurationFromAnchorWvs_Sg[1];
						DurationFromAnchorWvs_Sg[1]=DurationFromAnchorWvs_Sg[1]+Offset_segment[1];

						for (m=NumberOfSegmentsWws[j];m>1; m--) 
						{
							Offset_segment[m]=LastSample-(int)NumberOfSamplesWvs_Sg[m]-DurationFromAnchorWvs_Sg[m];
							DurationFromAnchorWvs_Sg[m]=DurationFromAnchorWvs_Sg[m]+Offset_segment[m];
							LastSample=Offset_segment[m];
						}
					}	
					if (NumberOfSamplesWvs_Sg[NumberOfSegmentsWws[j]]<20) //if last segment has less than 20 samples
					{
						NumberOfSegmentsWws[j]=NumberOfSegmentsWws[j]-1; //reduce the number of segments output
					}
					if (NumberOfSamplesWvs_Sg[NumberOfSegmentsWws[j]]<20) //if the second last segment has less than 20 samples
					{
						NumberOfSegmentsWws[j]=NumberOfSegmentsWws[j]-1;
					}
				}
				if ((NumberOfSegmentsWws[j]==1)&&(Type[PulseDescrIndex+PulseDescrlocation][j]==2))// for received, check if there is offset between Duration From Achor and FirstSample
				{
					if (DurationFromAnchorWvs_Sg[1]<FirstSample)
					{
						Offset_segment[1]=FirstSample-DurationFromAnchorWvs_Sg[1];
						DurationFromAnchorWvs_Sg[1]=DurationFromAnchorWvs_Sg[1]+Offset_segment[1];
					}
				}

				for (m=1;m<=NumberOfSegmentsWws[j]; m++) 
				{
					if ((Type[PulseDescrIndex+PulseDescrlocation][j]!=1)) // for received waveforms only 
					{
				
						/* Calculate peak */
           				max_waveform_elevation=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]);
            			peak_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]);
	            		peak_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]);
   		    			peak_h=max_waveform_elevation;

						output_waveform=0;
						if (     
							( (area_type==1) &&
							(sqrt(
							(peak_e-easting_centre)*(peak_e-easting_centre)+
							(peak_n-northing_centre)*(peak_n-northing_centre))
							<radius))
							||
    
							/* If peak within square */ ((area_type==2)&&(peak_e>west_edge)&&(peak_e<east_edge)&&(peak_n>south_edge)&&(peak_n<north_edge))
							)
							output_waveform=1;

						if ((area_type==2)&&
							(peak_e>west_edge)&&(peak_e<east_edge)&&
							(peak_n>south_edge)&&(peak_n<north_edge))
							output_waveform=1;

						if (area_type==3)
						{
							for (k=1;k<=NumberOfSamplesWvs[j];k++)
							{
								sample_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+k);
	    	    				sample_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+k);
		    					sample_h=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+k);
								if((sample_e>west_edge)&&(sample_e<east_edge)
									&&(sample_n>south_edge)&&(sample_n<north_edge)
									&&(sample_h>low_edge)&&(sample_h<high_edge)
									&&(surface_return_waveform[k]>extract_threshold))
									output_waveform=1;

							}
						}
						if (area_type==0 || process_all==1)
							output_waveform=1;

						// For Riegl Q680 Check if there are two channels in the pulse descriptor for a given waveform
						// if high power channel has 8 bits only low power channel is exported

						if ( (NumberOfSamplings[PulseDescrIndex+PulseDescrlocation]==3)
							&&(Type[PulseDescrIndex+PulseDescrlocation][2]==2)
							&&(Type[PulseDescrIndex+PulseDescrlocation][3]==2)
							&&(Channel[PulseDescrIndex+PulseDescrlocation][2]==1) /*low power channel*/
							&&(Channel[PulseDescrIndex+PulseDescrlocation][3]==0))/*high power channel*/
						{
							if (BitsPerSample[PulseDescrIndex+PulseDescrlocation][3]==8)
							{
								if (j==2)
									output_waveform=1;
								if (j==3)
									output_waveform=0;
							}
							else if (BitsPerSample[PulseDescrIndex+PulseDescrlocation][3]==16)
							{
								if (j==2)
								{
									if (max_surface_return<255)
									{
									output_waveform=1;
									low_saturated=0; 
									}
									else
									{
									output_waveform=0;
									low_saturated=1; //low channel saturated
									}
								}
								if (j==3)
								{
									if (low_saturated=1)
										output_waveform=1;
									else 
										output_waveform=0;
								}
							}
						}

						////////////////////////////////////////////////////
						//CHECK for anomalies if within the area of interest
						if (output_waveform==1)
						{
							// For Riegl waveforms that have less than 20 bins of returning waveform and it is the only segment - to eliminate them
							// or if the offset to waveorm gives a negative elevation
							if ((start_pulse_elevation+waveform_vector_elevation*(int)DurationFromAnchorWvs_Sg[m]<0&&NumberOfSegmentsWws[j]==1)
								||(NumberOfSegmentsWws[j]==1&&NumberOfSamplesWvs_Sg[m]<=20))
							{
								output_waveform=0;
							}
							//if there are two segments but the first one is too short (<20) then disregard 
							if (NumberOfSegmentsWws[j]>1)
							{ 
								if (m==1)
								{
									if( NumberOfSamplesWvs_Sg[m]<=20)
									{
										output_waveform=0;
										segment_too_short=1;
									}
									else
									{
										segment_too_short=0;
									}
								}
								if (m>1&&segment_too_short==1)
								{
									output_waveform=0;
								}
							}
							//if a waveform has more than one sampling and the first sampling is the trasmitted waveform
							//if that transmitted waveform has less than 30 samples the whole waveform will be disregarded
							if (NumberOfSamplings[PulseDescrIndex+PulseDescrlocation]>1)
							{
								if ((Type[PulseDescrIndex+PulseDescrlocation][1]==1) 
									&&((NumberOfSamplesWvs_Sg[1]<20)
									||(max_start_pulse_bin_Sg[1]<=5)
									||(max_start_pulse_bin_Sg[1]==NumberOfSamplesWvs_Sg[1])
									||(max_start_pulse_bin_Sg[1]==NumberOfSamplesWvs_Sg[1]-1))) //assuming there is only one transmitted waveform segment
								{
									output_waveform=0;		  
								}
							}
							if ((output_waveform==0)&&(NumberOfSegmentsWws[j]==m))
								count_disregarded++;
						}

						if (output_waveform==1)
						{

    						/* Find and output waveform */
								  
							if (Type[PulseDescrIndex+PulseDescrlocation][j]==2&&m==1) //if there is returning waveform and it is its first segment, output the outgoing waveform information (to avoid repetition for the second segment)
							{
								//outgoing waveform
								previous_gps_time=GPSTime;
	    		    			fprintf(out_waveform_fp,"%d,%.0lf,%d,%d,%lf,%lf,%lf,",i,previous_gps_time,(int)SampleUnitsComposition[PulseDescrIndex+PulseDescrlocation],PulseDescrIndex,start_pulse_easting,start_pulse_northing,start_pulse_elevation);

								no_transmitted=0;
								if (Type[PulseDescrIndex+PulseDescrlocation][j]==2&&NumberOfSamplings[PulseDescrIndex+PulseDescrlocation]==1)
								{
									fprintf(out_waveform_fp,"0,"); // no transmitted waveform recorded
									no_transmitted=1;
								}
								else if (start_pulse_length_Sg[m]>=98)
								{
									fprintf(out_waveform_fp,"98,");
									for (k=1;k<=98;k++)
										fprintf(out_waveform_fp,"%d,",start_pulse_waveform_Sg[m][k]);
								}
								else
								{
									fprintf(out_waveform_fp,"%d,",start_pulse_length_Sg[m]);
									for (k=1;k<=start_pulse_length_Sg[m];k++)
										fprintf(out_waveform_fp,"%d,",start_pulse_waveform_Sg[m][k]);
								}

								fprintf(out_waveform_fp,"%lf,%lf,%lf,",waveform_vector_easting, waveform_vector_northing, waveform_vector_elevation);
								fprintf(out_waveform_fp,"%d,",NumberOfSegmentsWws[j]);
							}

							fprintf(out_waveform_fp,"%d,", m); //segment number

							// E,N,E coordinates of the begining and end of the waveform (begining of the sampling)
							waveform_start_easting=start_pulse_easting+waveform_vector_easting*(int)DurationFromAnchorWvs_Sg[m];
    						waveform_start_northing=start_pulse_northing+waveform_vector_northing*(int)DurationFromAnchorWvs_Sg[m];
	    					waveform_start_elevation=start_pulse_elevation+waveform_vector_elevation*(int)DurationFromAnchorWvs_Sg[m];

							waveform_end_easting=waveform_start_easting+waveform_vector_easting*LastSample1;
							waveform_end_northing=waveform_start_northing+waveform_vector_northing*LastSample1;
							waveform_end_elevation=waveform_start_elevation+waveform_vector_elevation*LastSample1;
						  

							if (NumberOfSamplesWvs_Sg[m]<bin_truncation_limit)
								fprintf(out_waveform_fp,"%d,",(int)DurationFromAnchorWvs_Sg[m]);
							else
							{
								if (max_surface_bin>NumberOfSamplesWvs_Sg[m]-bin_truncation_limit)
									fprintf(out_waveform_fp,"%d,",(int)DurationFromAnchorWvs_Sg[m]+NumberOfSamplesWvs_Sg[m]-bin_truncation_limit);
								else
								{
									if (max_surface_bin<bin_truncation_limit)
										fprintf(out_waveform_fp,"%d,",(int)DurationFromAnchorWvs_Sg[m]);
									else
										fprintf(out_waveform_fp,"%d,",(int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]-bin_truncation_limit/2+1);
								}
							}

							fprintf(out_waveform_fp,"%lf,%lf,%lf,",waveform_start_easting, waveform_start_northing, waveform_start_elevation);
								  
							//output returning waveform segment

							/* Kludgy little bit to only output first 98 bins, to avoid overlong waveforms in excel */
	        				/* Revised to output number of bins, then series of bins, up to bin_truncation_limit, about 300 */
							if (NumberOfSamplesWvs_Sg[m]<bin_truncation_limit)
							{
								fprintf(out_waveform_fp,"%d,",NumberOfSamplesWvs_Sg[m]);
		    					for (k=1;k<=NumberOfSamplesWvs_Sg[m];k++)
    			    				fprintf(out_waveform_fp,"%d,",surface_return_waveform_Sg[m][k]);

    							for (k=1;k<=NumberOfSamplesWvs_Sg[m];k++)
	    						{
									sample_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+k);
	    	    					sample_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+k);
		    						sample_h=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+k);
    							}
							}
							else
							{
								if (max_surface_bin_Sg[m]>NumberOfSamplesWvs_Sg[m]-bin_truncation_limit)
								{
									fprintf(out_waveform_fp,"%d,",bin_truncation_limit);
		    						for (k=NumberOfSamplesWvs_Sg[m]-bin_truncation_limit+1;k<=NumberOfSamplesWvs_Sg[m];k++)
    			    					fprintf(out_waveform_fp,"%d,",surface_return_waveform_Sg[m][k]);

    								for (k=NumberOfSamplesWvs_Sg[m]-bin_truncation_limit+1;k<=NumberOfSamplesWvs_Sg[m];k++)
	    							{
										sample_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+k);
	    	    						sample_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+k);
		    							sample_h=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+k);
    								}
								}
								else
								{
									if (max_surface_bin_Sg[m]<bin_truncation_limit)
									{
										/* output the start */
										fprintf(out_waveform_fp,"%d,",bin_truncation_limit);
		    							for (k=1;k<=bin_truncation_limit;k++)
    			    						fprintf(out_waveform_fp,"%d,",surface_return_waveform_Sg[m][k]);

										for (k=1;k<=bin_truncation_limit;k++)
	    								{
											sample_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+k);
	    	    							sample_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+k);
		    								sample_h=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+k);
										}
									}
									else
									{
										/* Output around the peak */
										fprintf(out_waveform_fp,"%d,",bin_truncation_limit);

		    							for (k=max_surface_bin_Sg[m]-bin_truncation_limit/2+1;k<=max_surface_bin_Sg[m]+bin_truncation_limit/2;k++)
    			    						fprintf(out_waveform_fp,"%d,",surface_return_waveform_Sg[m][k]);

										for (k=max_surface_bin_Sg[m]-bin_truncation_limit/2+1;k<=max_surface_bin_Sg[m]+bin_truncation_limit/2;k++)
	    								{
											sample_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+k);
	    	    							sample_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+k);
		    								sample_h=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+k);
										}
									}
								}
							}


            				max_waveform_elevation=start_pulse_elevation+waveform_vector_elevation*((int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]);
	            			peak_e=start_pulse_easting+waveform_vector_easting*((int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]);
		            		peak_n=start_pulse_northing+waveform_vector_northing*((int)DurationFromAnchorWvs_Sg[m]+max_surface_bin_Sg[m]);
    		    			peak_h=max_waveform_elevation;
			  
							if (m==NumberOfSegmentsWws[j])
							{
							fprintf(out_waveform_fp,"%lf\n",max_waveform_elevation); 
							waveform_count++;
							}
							else
								fprintf(out_waveform_fp,"%lf,",max_waveform_elevation); 

    						if (NumberOfSamplesWvs_Sg[m]<=5000)
	    						waveform_bin_stats[(int)NumberOfSamplesWvs_Sg[m]]++;
		        					

							/* Check flightline extremes */
							if ((start_pulse_easting<westmost_aircraft_e)||(westmost_aircraft_e<-900))
							{
								westmost_aircraft_e=start_pulse_easting;
								westmost_aircraft_n=start_pulse_northing;
							}

							if ((start_pulse_easting>eastmost_aircraft_e)||(eastmost_aircraft_e<-900))
							{
								eastmost_aircraft_e=start_pulse_easting;
								eastmost_aircraft_n=start_pulse_northing;
							}
						} //end of the waveforms within specified area
					} // output only for received waveforms (includes outgoing)
					if (NumberOfSamplings[PulseDescrIndex+PulseDescrlocation]==1&&Type[PulseDescrIndex+PulseDescrlocation][j]==1)
					{
						count_transmittedonly++; //count pulses with no returning waveform
					}
					if (NumberOfSegmentsWws[j]>=2&&m==NumberOfSegmentsWws[j])
						count_multisegment++;
				}// End of number of segments loop
			}// End of number of samplings loop (Sampling 1=outgoing; Sampling 2=returning; if more than 2 it will increase waveform shot count) j-loop	
			} // End of pulse records loop (number of laser shots) i-loop

			fprintf(log_fp,"\n* Waveform processed: \t\t\t\t%d *\n", NumberOfPulses);
			fprintf(log_fp,"* Number of multisegment waveforms: \t\t%d *\n", count_multisegment);
			fprintf(log_fp,"* Waveforms with transmitted pulse only: \t%d *\n", count_transmittedonly);
				  
			fprintf(log_fp,"* Waveforms in selected area: \t\t\t%d *\n", waveform_count);
					fprintf(log_fp,"\n\tNumber of recorded samples \tNumber of waveforms\n");
			  
			for (i=1;i<=5000;i++)
				if (waveform_bin_stats[i]>0)
					fprintf(log_fp,"\t\t%d \t\t\t\t%d\n",i,waveform_bin_stats[i]);		

			fprintf(log_fp,"\n\t Writing to file: <<%s>>\n",out_waveform_filename);

			fclose(out_waveform_fp);

		} /* IF - check if point within extent */
		else
		{
			printf("\n%s - None\n",filename);
			fprintf(log_fp,"\nNo waveforms within specified area: %s - None\n",filename);
		}
			  
	if (no_transmitted==1)
		fprintf(log_fp,"\nWARNING! Some or all waveforms in the data do not contain recording of \nthe transmitted pulse!\n");
	
	fprintf(log_fp,"\n***********************************************************************");
	fprintf(log_fp,"\n\nPROCESSING COMPLETED!\n");
	printf("PROCESSING COMPLETE!\n");

	return(0);
}