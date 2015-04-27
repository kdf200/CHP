//***********************************************************************
//**               Module 3 VIS of Canopy Height Program               **
//***********************************************************************
//** Isolates ground points from single and last returns.              **
//** Provides visualization.                                            **
//**                                                                   **
//** Authors: Karolina Fieber, Ian Davenport                           **
//** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **
//** University of Reading, UK              November 2014 - March 2015 **
//***********************************************************************

#include<math.h>
#include<string.h>
#include<stdio.h>
#include<malloc.h>
#include "plcdemos.h"

static int ns      = 20; 
static int colorbar = 1,colorbar1 = 1;

void save_plot( char * );

double readindouble(FILE *fp)
{
	double scanned_number;
	char instring[200];
	int i;
	instring[0]=' ';
	for (i=1;(instring[i-1]!=',')&&(instring[i-1]!='\n')&&(i<99);i++)
		instring[i]=fgetc(fp);
	if (i<99)
	{
		instring[i-1]='\0';
		sscanf(instring,"%lf",&scanned_number);
	}
	else 
		scanned_number=-999;
	return scanned_number;
}

int readininteger(FILE *fp)
{
	int scanned_number;
	char instring[200];
	int i;
	instring[0]=' ';
	for (i=1;(instring[i-1]!=',')&&(instring[i-1]!='\n')&&(i<99);i++)
		instring[i]=fgetc(fp);
	if (i<99)
	{
		instring[i-1]='\0';
		sscanf(instring,"%ld",&scanned_number);
	}
	else 
		scanned_number=-999;
	return scanned_number;
}

/************************************************************************************************************************************/
/*													P R O G R A M																	*/
/************************************************************************************************************************************/
int
main( int argc, const char *argv[] )
{
	/*****************/
	/* File pointers */
	/*****************/

	FILE *in_fp1,*in_fp2;
	FILE *out_fp1,*out_fp2, *in_len;
	FILE *out_sum;
	
	/*************/
	/* Variables */
	/*************/

	char in_filename[200],in_filename2[200];
    char out_filename1[200],out_filename2[200];
	char out_filename_sum[200],in_filename_len[200];

	char	f_name_single[25]="_ground_single.psc"; 
	char	f_name_last[25]="_ground_last.psc"; 
	char	f_name_single1[250], f_name_last1[250];
	
	int count_points;
	int total_number_of_points;
	int waveform_one_peak;
	int larger_radius;

	/* New variables for 2011 version of data */


	//allowing ten segments in waveform sampling

	char identifier[25];
	int Sampling_unit;
	double mean_width_transmitted,noise_estim;
	double mean_vector_elevation;

	// DTM reference
	int  i,j,k,no_good, count_DTM, A;
	int *singleID, *singleNo, *singleAmplitude;
	double  mean_single_amplitude,mean_single_width;
	double E,N,Z, ground_E, ground_N, ground_Z, cum_Z_DTM, mean_Z_DTM,distance_DTM, W,min_Z;
	double *singleWidth, *singleEasting, *singleNorthing, *singleElevation;

	double mean_amplitude_trans;
	double sum_Z_DTM_wgt,sum_weight,mean_Z_DTM_wgt,distance,weight;
	int steep,length_data;


	//DTM last
	int *lastID, *lastNo, *lastAmplitude, *lastNoAll;
	double *lastWidth, *lastEasting, *lastNorthing, *lastElevation, *lastAngle;

	double E_min, E_max, N_min, N_max, Z_min, Z_max;

	/// colorbar
	double      *clevel,*clevel1, *shedge, *shedge1;
	double      fill_width = 2., cont_width = 0.;
    double      colorbar_width, colorbar_height;
    int      cont_color = 0;
#define NUM_AXES    1
    int      n_axis_opts  = NUM_AXES;
    const char *axis_opts[] = {
        "bcvtm",
    };
    int      num_values[NUM_AXES],num_values1[NUM_AXES];
    double      *values[NUM_AXES],*values1[NUM_AXES];
    double      axis_ticks[NUM_AXES] = {
        0.0,
    };
    int      axis_subticks[NUM_AXES] = {
        0,
    };
#define NUM_LABELS    1
    int      n_labels     = NUM_LABELS;
    int      label_opts[] = {
        PL_COLORBAR_LABEL_BOTTOM,
    };
    const char *labels[] = {
        "Elevation [m]",
    };

	char         string[30];


	E_min=-999; E_max=-999; N_min=-999; N_max=-999;Z_min=-999; Z_max=-999;


	/*****************/
	/* P R O G R A M */
	/*****************/

	printf("***********************************************************************\n");
	printf("**               Module 3 VIS of Canopy Height Program               **\n");
	printf("***********************************************************************\n");
	printf("** It also isolates ground points from single and last returns.      **\n");
	printf("** Provides visualization.                                           **\n");
	printf("**                                                                   **\n");
	printf("** Authors: Karolina Fieber, Ian Davenport                           **\n");
	printf("** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	printf("** University of Reading, UK                              March 2015 **\n");
	printf("***********************************************************************\n\n");

	printf("Enter the name of the input file (without '_fulwvs' suffix): ");
	scanf("%s",&identifier);

	strcpy(f_name_single1,identifier);
	strcpy(f_name_last1,identifier);

	strcat(f_name_single1,f_name_single);
	strcat(f_name_last1,f_name_last);
	
	sprintf(out_filename_sum,"%s_info_Module3.txt",identifier);
	out_sum=fopen(out_filename_sum,"w"); 
	if (!out_sum)
	{
		printf("Can't open output data file <%s>\n",out_filename_sum);
 		return(-1);
	}

	fprintf(out_sum,"***********************************************************************\n");
	fprintf(out_sum,"**               Module 3 VIS of Canopy Height Program               **\n");
	fprintf(out_sum,"***********************************************************************\n");
	fprintf(out_sum,"** It also isolates ground points from single and last returns.      **\n");
	fprintf(out_sum,"** Provides visualization.                                           **\n");
	fprintf(out_sum,"**                                                                   **\n");
	fprintf(out_sum,"** Authors: Karolina Fieber, Ian Davenport                           **\n");
	fprintf(out_sum,"** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	fprintf(out_sum,"** University of Reading, UK                              March 2015 **\n");
	fprintf(out_sum,"***********************************************************************\n\n");
	//fprintf(out_sum,"Processing file <<%s>>\n",in_filename);
	
	printf("\nComparison of: \n\t1- mean of point elevation (flat terrain) \n\t2- weighted average of point elevation (steep terrain)?  ");
	scanf("%d",&steep);
	while (steep!=1 && steep!=2)
	{
		printf("\n\tWrong input: 1 or 2: ");
		scanf("%lf",&steep);
	}
	fprintf(out_sum, "Comparison of: \n\t1- mean of point elevation (flat terrain) \n\t2- weighted average of point elevation (steep terrain)? %d",steep);

	printf("\nSpecify the search area for mean ground elevation estimation in [m]\n\t(in selection of ground points): +- ");
	scanf("%lf",&distance_DTM);
	while (distance_DTM<=0 && distance_DTM>=100)
	{
		printf("\n\tWrong input: must be larger than 0: ");
		scanf("%lf",&distance_DTM);
	}
	fprintf(out_sum, "\n\nArea used for ground point selection [m]: +- %.2lf\n",distance_DTM);
	fprintf(out_sum, "\n***********************************************************************\n");
	printf("\n***********************************************************************\n");

	// OPEN file with parameters
	sprintf(in_filename_len,"%s_param.txt",identifier);
	in_len=fopen(in_filename_len,"r");
	if (!in_len)
	{
		printf("\nCan't open input data file <%s>\n",in_filename_len);
		fprintf(out_sum,"\n\tCan't open input data file <%s>\n",in_filename_len);
		fclose(out_sum);
 		return(-1);
	}
	else
	{
		printf("\n\tReading from file <<%s>>\n",in_filename_len);
		fprintf(out_sum,"\n\tReading from file <<%s>>\n",in_filename_len);
	}

	// read in parameters
	length_data=readininteger(in_len);
	waveform_one_peak=readininteger(in_len);
	total_number_of_points=readininteger(in_len);
	Sampling_unit=readininteger(in_len);
	mean_amplitude_trans=readindouble(in_len);
	mean_width_transmitted=readindouble(in_len);
	mean_single_amplitude=readindouble(in_len);
	mean_single_width=readindouble(in_len);
	noise_estim=readindouble(in_len);
	mean_vector_elevation=readindouble(in_len);

	singleID=(int *)calloc(total_number_of_points+2,sizeof(int));
	singleNo=(int *)calloc(total_number_of_points+2,sizeof(int));
	singleAmplitude=(int *)calloc(total_number_of_points+2,sizeof(int));	
	singleWidth=(double *)calloc(total_number_of_points+2,sizeof(double));
	singleEasting=(double *)calloc(total_number_of_points+2,sizeof(double));
	singleNorthing=(double *)calloc(total_number_of_points+2,sizeof(double));
	singleElevation=(double *)calloc(total_number_of_points+2,sizeof(double));


	lastID=(int *)calloc(total_number_of_points+2,sizeof(int));
	lastNo=(int *)calloc(total_number_of_points+2,sizeof(int));
	lastNoAll=(int *)calloc(total_number_of_points+2,sizeof(int));
	lastAmplitude=(int *)calloc(total_number_of_points+2,sizeof(int));
	lastWidth=(double *)calloc(total_number_of_points+2,sizeof(double));
	lastEasting=(double *)calloc(total_number_of_points+2,sizeof(double));
	lastNorthing=(double *)calloc(total_number_of_points+2,sizeof(double));
	lastElevation=(double *)calloc(total_number_of_points+2,sizeof(double));
	lastAngle=(double *)calloc(total_number_of_points+2,sizeof(double));

	larger_radius=0;
	///////////////////////////////////////////////////////////////////
	// Ground point selection
	///////////////////////////////////////////////////////////////////

	//read in single 

	printf("\n***********************************************************************");
	printf("\n**********  Ground point selection from single returns  ***************");
	printf("\n***********************************************************************");
	          
	fprintf(out_sum,"\n***********************************************************************");
	fprintf(out_sum,"\n**********  Ground point selection from single returns  ***************");
	fprintf(out_sum,"\n***********************************************************************");                 

	sprintf(in_filename,"%s_XYZ_single.txt",identifier);
	in_fp1=fopen(in_filename,"r");
	if (!in_fp1)
	{
		printf("Can't open input data file <%s>\n",in_filename);
		fprintf(out_sum,"\n\tCan't open input data file <%s>\n",in_filename);
 		return(-1);
	}
	else
	{
		printf("\n\n\tReading from file <<%s>>\n",in_filename);
		fprintf(out_sum,"\n\n\tReading from file <<%s>>\n",in_filename);
	}

	sprintf(out_filename1,"%s_ground_single.txt",identifier);
	out_fp1=fopen(out_filename1,"w");
	if (!out_fp1)
	{
		printf("Can't open output data file <%s>\n",out_filename1);
		fprintf(out_sum,"\n\tCan't open input data file <%s>\n",out_filename1);
 		return(-1);
	}
	else
	{
		printf("\tWriting to file <<%s>>\n\n",out_filename1);
		fprintf(out_sum,"\tWriting to file <<%s>>\n",out_filename1);
	}


	for (i=1;i<=waveform_one_peak;i++)
	{
		singleID[i]=readininteger(in_fp1); //pulseID
		singleNo[i]=readininteger(in_fp1); //pulse number
		singleAmplitude[i]=readininteger(in_fp1); //amplitude
		singleWidth[i]=readindouble(in_fp1); //width
		singleEasting[i]=readindouble(in_fp1); //easting  
		singleNorthing[i]=readindouble(in_fp1);//northing
		singleElevation[i]=readindouble(in_fp1); //elevation
	}

	///////////////////////////////////
	//// FIRST Smoothing iteration ////
	///////////////////////////////////

	for (i=1;i<=waveform_one_peak;i++) // for every point in the file check the points within its set proximity for their elevation values
	{
		E=singleEasting[i];
		N=singleNorthing[i];
		Z=singleElevation[i];
		W=singleWidth[i];
		A=singleAmplitude[i];
		no_good=0;
		cum_Z_DTM=0;
		count_DTM=0;
		min_Z=Z;

		sum_Z_DTM_wgt=0;
		sum_weight=0;
		mean_Z_DTM_wgt=0;

		for(j=1;j<=waveform_one_peak;j++)
		{
			ground_E=singleEasting[j];
			ground_N=singleNorthing[j];
			ground_Z=singleElevation[j];
			
			if (((ground_E<=(E+distance_DTM))&&(ground_E>=(E-distance_DTM)))
				&&((ground_N<=(N+distance_DTM))&&(ground_N>=(N-distance_DTM)))
				&&(i!=j)) 
			{
				if (steep==2)
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					distance=sqrt(((ground_E-E)*(ground_E-E))+((ground_N-N)*(ground_N-N)));
					weight=1/(distance*distance);
					sum_Z_DTM_wgt=sum_Z_DTM_wgt+(weight*ground_Z);
					sum_weight=sum_weight+weight;

					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
					}
				}
				else 
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
					}
				}
			}
		}

		mean_Z_DTM=cum_Z_DTM/count_DTM;
		mean_Z_DTM_wgt=sum_Z_DTM_wgt/sum_weight;

		if (count_DTM<1)
		{
			printf("Not enough points: specify larger radius!\n");
			larger_radius=1;
		}
		
		if (steep==2)
		{
			if (((Z<=(mean_Z_DTM_wgt+0.5))&&(A>mean_single_amplitude-mean_single_amplitude/4) && (W<=mean_width_transmitted+mean_width_transmitted/10))
				||(Z<mean_Z_DTM_wgt+0.2))
			{
				Z=Z;
			}
			else if (count_DTM>=4&&(Z>(mean_Z_DTM_wgt+0.5))
			||((Z>mean_Z_DTM_wgt+0.2)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
					if (((mean_Z_DTM_wgt-min_Z)>1)&&((Z-min_Z)>1.2))
					{
						Z=min_Z+(mean_Z_DTM_wgt-min_Z)/2; 
						singleElevation[i]=min_Z+(mean_Z_DTM_wgt-min_Z)/2;
					}
					else
					{
						Z=mean_Z_DTM_wgt; 
						singleElevation[i]=mean_Z_DTM_wgt;
					}
			}
		}
		else
		{
			if (count_DTM>=1&&(Z>(mean_Z_DTM+0.2))
			||((Z>mean_Z_DTM)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				if ((mean_Z_DTM-min_Z)>0.5)
				{
					Z=min_Z+(mean_Z_DTM-min_Z)/2; 
					singleElevation[i]=min_Z+(mean_Z_DTM-min_Z)/2;
				}
				else
				{
					Z=mean_Z_DTM; 
					singleElevation[i]=mean_Z_DTM;
				}
			}
		}
	}

	////////////////////////////////////
	//// SECOND Smoothing iteration ////
	////////////////////////////////////

	for (i=1;i<=waveform_one_peak;i++)
	{
		E=singleEasting[i];
		N=singleNorthing[i];
		Z=singleElevation[i];
		W=singleWidth[i];
		A=singleAmplitude[i];
		no_good=0;
		cum_Z_DTM=0;
		count_DTM=0;
		min_Z=Z;

		sum_Z_DTM_wgt=0;
		sum_weight=0;
		mean_Z_DTM_wgt=0;

		for(j=1;j<=waveform_one_peak;j++)
		{
			ground_E=singleEasting[j];
			ground_N=singleNorthing[j];
			ground_Z=singleElevation[j];
			
			if ((ground_E<=(E+distance_DTM)&&ground_E>=(E-distance_DTM))
				&&(ground_N<=(N+distance_DTM)&&ground_N>=(N-distance_DTM))
				&&(i!=j)) 
			{
				if (steep==2)
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					distance=sqrt(((ground_E-E)*(ground_E-E))+((ground_N-N)*(ground_N-N)));
					weight=1/(distance*distance);
					sum_Z_DTM_wgt=sum_Z_DTM_wgt+(weight*ground_Z);
					sum_weight=sum_weight+weight;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
					}
				}
				else
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
					}
				}
			}
		}

		mean_Z_DTM=cum_Z_DTM/count_DTM;
		mean_Z_DTM_wgt=sum_Z_DTM_wgt/sum_weight;


		if (steep==2)
		{
			if (((Z<=(mean_Z_DTM_wgt+0.2))&&(A>mean_single_amplitude-mean_single_amplitude/4) && (W<=mean_width_transmitted+mean_width_transmitted/10))
				||(Z<mean_Z_DTM_wgt+0.1))
			{
				Z=Z;
			}	
			else if (count_DTM>=4&&(Z>(mean_Z_DTM_wgt+0.2))
				||((Z>mean_Z_DTM_wgt+0.1)&&
				((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				{
					if (((mean_Z_DTM_wgt-min_Z)>0.8)&&((Z-min_Z)>0.8))
					{
						Z=min_Z+(mean_Z_DTM_wgt-min_Z)/2; 
						singleElevation[i]=min_Z+(mean_Z_DTM_wgt-min_Z)/2;
					}
					else
					{
						Z=mean_Z_DTM_wgt; 
						singleElevation[i]=mean_Z_DTM_wgt;
					}
				}
			}
		}
		else
		{
			if (count_DTM>=1&&(Z>(mean_Z_DTM+0.1))
			||((Z>mean_Z_DTM)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{

				if ((mean_Z_DTM-min_Z)>0.5)
				{
					Z=min_Z+(mean_Z_DTM-min_Z)/2; 
					singleElevation[i]=min_Z+(mean_Z_DTM-min_Z)/2;
				}
				else
				{
					Z=mean_Z_DTM; 
					singleElevation[i]=mean_Z_DTM;
				}
			}
		}
	}

	///////////////////////////////////
	//// THIRD Smoothing iteration ////
	///////////////////////////////////

	for (i=1;i<=waveform_one_peak;i++)
	{
		E=singleEasting[i];
		N=singleNorthing[i];
		Z=singleElevation[i];
		W=singleWidth[i];
		A=singleAmplitude[i];
		no_good=0;
		cum_Z_DTM=0;
		count_DTM=0;

		sum_Z_DTM_wgt=0;
		sum_weight=0;
		mean_Z_DTM_wgt=0;

		//find the minimum and maximum
		if ((E<E_min)||(E_min==-999))
			E_min=E;

		if ((E>E_max)||(E_max==-999))
			E_max=E;

		if ((N<N_min)||(N_min==-999))
			N_min=N;

		if ((N>N_max)||(N_max==-999))
			N_max=N;



		for(j=1;j<=waveform_one_peak;j++)
		{
			ground_E=singleEasting[j];
			ground_N=singleNorthing[j];
			ground_Z=singleElevation[j];
			
			if ((ground_E<=(E+distance_DTM)&&ground_E>=(E-distance_DTM))
				&&(ground_N<=(N+distance_DTM)&&ground_N>=(N-distance_DTM))
				&&(i!=j)) 
			{
				if (steep==2)
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					distance=sqrt(((ground_E-E)*(ground_E-E))+((ground_N-N)*(ground_N-N)));
					weight=1/(distance*distance);
					sum_Z_DTM_wgt=sum_Z_DTM_wgt+(weight*ground_Z);
					sum_weight=sum_weight+weight;
				}
				else
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
				}
			}
		}

		mean_Z_DTM=cum_Z_DTM/count_DTM;
		mean_Z_DTM_wgt=sum_Z_DTM_wgt/sum_weight;
		

		if (steep==2)
		{
			if (((Z<=(mean_Z_DTM_wgt+0.1))&&(A>mean_single_amplitude-mean_single_amplitude/4) && (W<=mean_width_transmitted+mean_width_transmitted/10))
				||(Z<mean_Z_DTM_wgt+0.05))
			{
				Z=Z;
			}
			else if (count_DTM>=4&&(Z>(mean_Z_DTM_wgt+0.1))
				||((Z>mean_Z_DTM_wgt+0.05)&&
				((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				Z=mean_Z_DTM_wgt;
				singleElevation[i]=mean_Z_DTM_wgt;
			}
		}
		else //flat
		{
			if (count_DTM>=1&&(Z>(mean_Z_DTM+0.05))
			||((Z>mean_Z_DTM)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				Z=mean_Z_DTM; 
				singleElevation[i]=mean_Z_DTM;
			}
		}


		if ((Z<Z_min)||(Z_min==-999))
			Z_min=Z;

		if ((Z>Z_max)||(Z_max==-999))
			Z_max=Z;

		fprintf(out_fp1,"%lf,",E);
		fprintf(out_fp1,"%lf,",N);
		fprintf(out_fp1,"%lf\n",Z);
	}

	printf("\n\tZmin=%.3lf, Z_max=%.3lf, dif= %.3lf\n", Z_min, Z_max, Z_max-Z_min);
	
	fclose(in_fp1);
	fclose(out_fp1);
	fprintf(in_len,"%d\n",waveform_one_peak);
	fprintf(out_sum,"\tGround points from single returns: %d\n",waveform_one_peak);
	 if (larger_radius==1)
	 {
		 fprintf(out_sum,"\n\tWARNING: larger radius might be needed!\n");
	 }
	 larger_radius=0;


	//////////////////////////////////////////////////////////////////////
	/// VISUALIZATION
	/////////////////////////////////////////////
	plspal0( "cmap0_black_on_white.pal" );
	plsdev( "wingcc" );
	plinit();
	plssub( 1, 2);

	pladv( 0 );
    plcol0( 1 );
    plvpor( 0.0, 1.0, 0.0, 0.9 );
    plwind( -1.0, 1.0, -1.0, 1.0 );
	plschr( 0., 0.7 );
    plw3d( 1.0, ((N_max-N_min)/(E_max-E_min)), ((Z_max-Z_min)/(E_max-E_min)), E_min, E_max, N_min, N_max, Z_min, Z_max, 50.0, 45. );
    plbox3( "bnstu", "Easting [m]", 0.0, 0,
        "bnstu", "Northing [m]", 0.0, 0,
        "bcdtu", "Elevation [m]", 0.0, 0 );


	plspal0( "cmap0_black_on_white.pal" );
	plspal1( "cmap1_blue_red.pal", 1 );

	for ( i = 1; i <= waveform_one_peak; i++ )
    {
		plcol1( ( singleElevation[i] - Z_min ) / ( Z_max - Z_min ) );
		plstring3( 1, &singleEasting[i], &singleNorthing[i],&singleElevation[i],"#(729)"  );
	 }
	plschr( 0., 1.0 );
	plcol0( 1 );
	plmtex( "t", 0.1, 0.5, 0.5, "SINGLE: (alt: 50; az: 45)" );

	////////////////////////////////////////////////////////////////////
	// Create a labelled box to hold the plot.
	clevel = (double *) calloc( (size_t) ns, sizeof ( double ) );
	shedge = (double *) calloc( (size_t) ( ns + 1 ), sizeof ( double ) );

	for ( i = 0; i < ns; i++ )
		clevel[i] = Z_min + ( Z_max- Z_min) * ( i + 0.5 ) / (double) ns;

	for ( i = 0; i < ns + 1; i++ )
		shedge[i] = Z_min + ( Z_max - Z_min ) * (double) i / (double) ns;

	plschr( 0., 0.5 );
	plspal0( "cmap0_black_on_white.pal" );
	plspal1( "cmap1_blue_red.pal", 1 );

	plcol0( 1 );
	plschr( 0., 0.5 );
    plenv( E_min, E_max, N_min, N_max, 1, 0 );

    // Plot the data that was prepared above.
	 for ( i = 1; i <= waveform_one_peak; i++ )
    {
		plcol1( ( singleElevation[i] - Z_min ) / ( Z_max - Z_min ) );
		plstring( 1, &singleEasting[i], &singleNorthing[i],"#(729)"  );
	 }
	plschr( 0., 0.5 );

	// draw colorbar
	if ( colorbar )
	{
		// Smaller text
		plschr( 0.0, 0.75 );
		// Small ticks on the vertical axis
		plsmaj( 0.0, 0.5 );
		plsmin( 0.0, 0.5 );

		//printf("\nGenerating colorbar");
		num_values[0] = ns + 1;
		values[0]     = shedge;
		plcolorbar( &colorbar_width, &colorbar_height,
			PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL, 0,
			0.005, 0.0, 0.0375, 0.875, 0, 1, 1, 0.0, 0.0,
			cont_color, cont_width,
			n_labels, label_opts, labels,
			n_axis_opts, axis_opts,
			axis_ticks, axis_subticks,
			num_values, (const double * const *) values );

		// Reset text and tick sizes
		plschr( 0.0, 1.0 );
		plsmaj( 0.0, 1.0 );
		plsmin( 0.0, 1.0 );
	}
	plschr( 0.,1.0 );
	plmtex("b",1.7, 0.5, 0.5, "Easting [m]");
	plmtex("l",3.5, 0.5, 0.5, "Northing [m]");
	plmtex("t",1.5, 0.5, 0.5, "SINGLE: Plan view ");
	save_plot( f_name_single1 );

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	////DTM from last
	//read in last
	printf("\n***********************************************************************");
	printf("\n************  Ground point selection from last returns  ***************");
	printf("\n***********************************************************************");

	fprintf(out_sum,"\n***********************************************************************");
	fprintf(out_sum,"\n*************  Ground point selection from last returns  **************");
	fprintf(out_sum,"\n***********************************************************************");

	E_min=-999; E_max=-999; N_min=-999; N_max=-999;Z_min=-999; Z_max=-999;


	sprintf(in_filename2,"%s_XYZ.txt",identifier);
	in_fp2=fopen(in_filename2,"r");
	if (!in_fp2)
	{
		printf("Can't open input data file <%s>\n",in_filename2);
		fprintf(out_sum,"\nCan't open input data file <%s>\n",in_filename2);
 		return(-1);
	}
	else
	{
		printf("\n\n\tReading from file <<%s>>\n",in_filename2);
		fprintf(out_sum,"\n\n\tReading from file <<%s>>\n",in_filename2);
	}

	sprintf(out_filename2,"%s_ground_last.txt",identifier);
	out_fp2=fopen(out_filename2,"w");
	if (!out_fp2)
	{
		printf("Can't open output data file <%s>\n",out_filename2);
 		return(-1);
	}
	else
	{
		printf("\tWriting to file <<%s>>\n\n",out_filename2);
		fprintf(out_sum,"\tWriting to file <<%s>>\n",out_filename2);
	}

	for (i=1;i<=total_number_of_points;i++)
	{
		lastID[i]=readininteger(in_fp2); //pulseID
		lastNo[i]=readininteger(in_fp2); //pulse number
		lastNoAll[i]=readininteger(in_fp2);
		lastAmplitude[i]=readininteger(in_fp2);//amplitude
		lastWidth[i]=readindouble(in_fp2);//width
		lastEasting[i]=readindouble(in_fp2); //easting
		lastNorthing[i]=readindouble(in_fp2);//northing
		lastElevation[i]=readindouble(in_fp2);//elevation
		lastAngle[i]=readindouble(in_fp2);
	}

	k=1;
	for (i=1;i<=total_number_of_points;i++)
	{
		if (lastNo[i]==lastNoAll[i]&&lastNo[i]==1) // for single peaks only
		{

			singleEasting[k]=lastEasting[i];
			singleNorthing[k]=lastNorthing[i];
			singleElevation[k]=lastElevation[i];
			singleWidth[k]=lastWidth[i];
			singleAmplitude[k]=lastAmplitude[i];
			k=k+1;
	
		}
		else if (lastNo[i]==lastNoAll[i]&&lastElevation[i]<lastElevation[i-1]) // for multi peaks when we want the last return (lowest elevation)
		{
			if (i==1)
			{
				singleEasting[k]=lastEasting[i];
				singleNorthing[k]=lastNorthing[i];
				singleElevation[k]=lastElevation[i];
				singleWidth[k]=lastWidth[i];
				singleAmplitude[k]=lastAmplitude[i];
				k=k+1;
			}
			else
			{
				singleEasting[k]=lastEasting[i];
				singleNorthing[k]=lastNorthing[i];
				singleElevation[k]=lastElevation[i];
				singleWidth[k]=lastWidth[i];
				singleAmplitude[k]=lastAmplitude[i];
				k=k+1;
			}
		}
		else if (lastNo[i]==lastNoAll[i]&&lastElevation[i]>lastElevation[i-1]) // for multi peaks when we want the first return (reverse recording)
		{
			if (i==1)
			{
				singleEasting[k]=lastEasting[i];
				singleNorthing[k]=lastNorthing[i];
				singleElevation[k]=lastElevation[i];
				singleWidth[k]=lastWidth[i];
				singleAmplitude[k]=lastAmplitude[i];
				k=k+1;
			}
			else
			{
				singleEasting[k]=lastEasting[i];
				singleNorthing[k]=lastNorthing[i];
				singleElevation[k]=lastElevation[i];
				singleWidth[k]=lastWidth[i];
				singleAmplitude[k]=lastAmplitude[i];
				k=k+1;
			}
		}
	}

	///////////////////////////////////
	//// FIRST Smoothing iteration ////
	///////////////////////////////////

	k=k-1;
	for (i=1;i<=k;i++) // for every point in the file check the points within its set proximity for their elevation values
	{
		E=singleEasting[i];
		N=singleNorthing[i];
		Z=singleElevation[i];
		W=singleWidth[i];
		A=singleAmplitude[i];
		no_good=0;
		cum_Z_DTM=0;
		count_DTM=0;
		min_Z=Z;

		sum_Z_DTM_wgt=0;
		sum_weight=0;
		mean_Z_DTM_wgt=0;

		for(j=1;j<=k;j++)
		{
			ground_E=singleEasting[j];
			ground_N=singleNorthing[j];
			ground_Z=singleElevation[j];
			
			if ((ground_E<=(E+distance_DTM)&&ground_E>=(E-distance_DTM))
				&&(ground_N<=(N+distance_DTM)&&ground_N>=(N-distance_DTM))
				&&(i!=j)) 
			{
				if (steep==2)
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					weight=1/(distance*distance);
					sum_Z_DTM_wgt=sum_Z_DTM_wgt+(weight*ground_Z);
					sum_weight=sum_weight+weight;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
	
					}
				}
				else 
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
	
					}
				}
			}
		}

		mean_Z_DTM=cum_Z_DTM/count_DTM;
		mean_Z_DTM_wgt=sum_Z_DTM_wgt/sum_weight;
		
		if (count_DTM<1)
		{
			printf("Not enough points: specify larger radius!\n");
			larger_radius=1;
		}

		if (steep==2)
		{
			if (count_DTM>=4&&(Z>(mean_Z_DTM_wgt+0.3))
				||((Z>mean_Z_DTM_wgt+0.3)&&
				((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				if ((mean_Z_DTM_wgt-min_Z)>0.5&&(mean_Z_DTM_wgt-min_Z)<1.0)
				{
					Z=min_Z+(mean_Z_DTM_wgt-min_Z)/2; 
					singleElevation[i]=min_Z+(mean_Z_DTM_wgt-min_Z)/2;
				}
				else if ((mean_Z_DTM_wgt-min_Z)>1.0)
				{
					Z=min_Z+(mean_Z_DTM_wgt-min_Z)/4; 
					singleElevation[i]=min_Z+(mean_Z_DTM_wgt-min_Z)/4;
				}
				else
				{
					Z=mean_Z_DTM_wgt; 
					singleElevation[i]=mean_Z_DTM_wgt;
				}
			}
			else if (((count_DTM>=4)&&(Z<(mean_Z_DTM_wgt-0.5))&&(W>mean_width_transmitted+mean_width_transmitted/8))&&
				(((A<3*noise_estim)&&((noise_estim/mean_amplitude_trans)<0.05))||((A<1.3*noise_estim)&&((noise_estim/mean_amplitude_trans)<0.12)))) // to avoid ringing echos
			{
				Z=mean_Z_DTM_wgt; 
				singleElevation[i]=mean_Z_DTM_wgt;
			}
		}
		else
		{
			if (count_DTM>=1&&(Z>(mean_Z_DTM+0.2))
			||((Z>mean_Z_DTM)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				if ((mean_Z_DTM-min_Z)>0.5&&(mean_Z_DTM-min_Z)<1.0)
				{
					Z=min_Z+(mean_Z_DTM-min_Z)/2; 
					singleElevation[i]=min_Z+(mean_Z_DTM-min_Z)/2;
				}
				else if ((mean_Z_DTM-min_Z)>1.0)
				{
					Z=min_Z+(mean_Z_DTM-min_Z)/4; 
					singleElevation[i]=min_Z+(mean_Z_DTM-min_Z)/4;
				}
				else
				{
					Z=mean_Z_DTM; 
					singleElevation[i]=mean_Z_DTM;
				}
			}
		}
	}

	////////////////////////////////////
	//// SECOND Smoothing iteration ////
	////////////////////////////////////

	for (i=1;i<=k;i++)
	{
		E=singleEasting[i];
		N=singleNorthing[i];
		Z=singleElevation[i];
		W=singleWidth[i];
		A=singleAmplitude[i];
		no_good=0;
		cum_Z_DTM=0;
		count_DTM=0;
		min_Z=Z;

		sum_Z_DTM_wgt=0;
		sum_weight=0;
		mean_Z_DTM_wgt=0;

		for(j=1;j<=k;j++)
		{
			ground_E=singleEasting[j];
			ground_N=singleNorthing[j];
			ground_Z=singleElevation[j];
			
			if ((ground_E<=(E+distance_DTM)&&ground_E>=(E-distance_DTM))
				&&(ground_N<=(N+distance_DTM)&&ground_N>=(N-distance_DTM))
				&&(i!=j)) 
			{
				if (steep==2)
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					weight=1/(distance*distance);
					sum_Z_DTM_wgt=sum_Z_DTM_wgt+(weight*ground_Z);
					sum_weight=sum_weight+weight;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
					}
				}
				else
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					if (min_Z>ground_Z)
					{
						min_Z=ground_Z;
					}
				}
			}
		}

		mean_Z_DTM=cum_Z_DTM/count_DTM;
		mean_Z_DTM_wgt=sum_Z_DTM_wgt/sum_weight;
		

		if (steep==2)
		{
			if (count_DTM>=4&&(Z>(mean_Z_DTM_wgt+0.2))
				||((Z>mean_Z_DTM_wgt+0.2)&&
				((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{

				if ((mean_Z_DTM_wgt-min_Z)>0.5)
				{
					Z=min_Z+(mean_Z_DTM_wgt-min_Z)/2; 
					singleElevation[i]=min_Z+(mean_Z_DTM_wgt-min_Z)/2;
				}
				else
				{
					Z=mean_Z_DTM_wgt; 
					singleElevation[i]=mean_Z_DTM_wgt;
				}
			}
			else if (((count_DTM>=4)&&(Z<(mean_Z_DTM_wgt-0.5))&&(W>mean_width_transmitted+mean_width_transmitted/8))&&
				(((A<3*noise_estim)&&((noise_estim/mean_amplitude_trans)<0.05))||((A<1.3*noise_estim)&&((noise_estim/mean_amplitude_trans)<0.12)))) // to avoid ringing echos
			{
				Z=mean_Z_DTM_wgt; 
				singleElevation[i]=mean_Z_DTM_wgt;
			}
		}
		else
		{
			if (count_DTM>=1&&(Z>(mean_Z_DTM+0.1))
			||((Z>mean_Z_DTM)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{

				if ((mean_Z_DTM-min_Z)>0.5)
				{
					Z=min_Z+(mean_Z_DTM-min_Z)/2; 
					singleElevation[i]=min_Z+(mean_Z_DTM-min_Z)/2;
				}
				else
				{
					Z=mean_Z_DTM; 
					singleElevation[i]=mean_Z_DTM;
				}
			}
		}
	}

	///////////////////////////////////
	//// THIRD Smoothing iteration ////
	///////////////////////////////////
	count_points=0;

	for (i=1;i<=k;i++)
	{
		E=singleEasting[i];
		N=singleNorthing[i];
		Z=singleElevation[i];
		W=singleWidth[i];
		A=singleAmplitude[i];
		no_good=0;
		cum_Z_DTM=0;
		count_DTM=0;

		sum_Z_DTM_wgt=0;
		sum_weight=0;
		mean_Z_DTM_wgt=0;

		//find the minimum and maximum
		if ((E<E_min)||(E_min==-999))
			E_min=E;

		if ((E>E_max)||(E_max==-999))
			E_max=E;

		if ((N<N_min)||(N_min==-999))
			N_min=N;

		if ((N>N_max)||(N_max==-999))
			N_max=N;


		for(j=1;j<=k;j++)
		{
			ground_E=singleEasting[j];
			ground_N=singleNorthing[j];
			ground_Z=singleElevation[j];
			
			if ((ground_E<=(E+distance_DTM)&&ground_E>=(E-distance_DTM))
				&&(ground_N<=(N+distance_DTM)&&ground_N>=(N-distance_DTM))
				&&(i!=j)) 
			{
				if (steep==2)
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
					weight=1/(distance*distance);
					sum_Z_DTM_wgt=sum_Z_DTM_wgt+(weight*ground_Z);
					sum_weight=sum_weight+weight;
				}
				else
				{
					cum_Z_DTM=cum_Z_DTM+ground_Z;
					count_DTM=count_DTM+1;
				}
			}
		}

		mean_Z_DTM=cum_Z_DTM/count_DTM;
		mean_Z_DTM_wgt=sum_Z_DTM_wgt/sum_weight;
		

		if (steep==2)
		{
			if (count_DTM>=4&&(Z>(mean_Z_DTM_wgt+0.1))
				||((Z>mean_Z_DTM_wgt+0.1)&&
				((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				Z=mean_Z_DTM_wgt; 
				singleElevation[i]=mean_Z_DTM_wgt;
			}
		}
		else
		{
			if (count_DTM>=1&&(Z>(mean_Z_DTM+0.05))
			||((Z>mean_Z_DTM)&&
			((A<mean_single_amplitude-mean_single_amplitude/3)||(W>mean_single_width)||(W>mean_width_transmitted+mean_width_transmitted/8))))
			{
				Z=mean_Z_DTM; 
				singleElevation[i]=mean_Z_DTM;
			}
		}

		if ((Z<Z_min)||(Z_min==-999))
			Z_min=Z;

		if ((Z>Z_max)||(Z_max==-999))
			Z_max=Z;

		fprintf(out_fp1,"%lf,",E);
		fprintf(out_fp1,"%lf,",N);
		fprintf(out_fp1,"%lf\n",Z);
	}
	
	printf("\n\tZmin= %.3lf, Zmax=%.3lf, dif=%.3lf\n", Z_min,Z_max,Z_max-Z_min);
	fprintf(out_sum, "\tGround points from last returns: %d\n",k);
	if (larger_radius==1)
	 {
		 fprintf(out_sum,"\n\tWARNING: larger radius might be needed!\n");
	 }
	fprintf(out_sum,"\n***********************************************************************");
	fprintf(out_sum, "\n\nPROCESSING COMPLETED!");
	printf( "\n\nPROCESSING COMPLETED!");
	
	fprintf(in_len,"%d\n",k);

	//////////////////////////////////////////////////////////////////////
	/// VISUALIZATION
	//////////////////////////////////////////////////////////////////////////////////////////
	pladv( 0 );
    plcol0( 1 );
    plvpor( 0.0, 1.0, 0.0, 0.9 );
    plwind( -1.0, 1.0, -1.0, 1.0 );
	plschr( 0.,0.7 );
    plw3d( 1.0, ((N_max-N_min)/(E_max-E_min)), ((Z_max-Z_min)/(E_max-E_min)), E_min, E_max, N_min, N_max, Z_min, Z_max, 50.0, 45. );
    plbox3( "bnstu", "Easting [m]", 0.0, 0,
        "bnstu", "Northing [m]", 0.0, 0,
        "bcdtu", "Elevation [m]", 0.0, 0 );

	for ( i = 1; i <= k-1; i++ )
    {
		plcol1( ( singleElevation[i] - Z_min ) / ( Z_max - Z_min ) );
		plstring3( 1, &singleEasting[i], &singleNorthing[i],&singleElevation[i],"#(729)" );
	 }
	plschr( 0.,1.0 );
	plcol0( 1 );
	plmtex( "t", 1.0, 0.5, 0.5, "LAST: (alt: 50; az: 45)" );

	//////////////////////////////////////////////////////////////////////////////////////////
	// Create a labelled box to hold the plot.
	clevel1 = (double *) calloc( (size_t) ns, sizeof ( double ) );
	shedge1 = (double *) calloc( (size_t) ( ns + 1 ), sizeof ( double ) );

	for ( i = 0; i < ns; i++ )
		clevel1[i] = Z_min + ( Z_max- Z_min) * ( i + 0.5 ) / (double) ns;

	for ( i = 0; i < ns + 1; i++ )
		shedge1[i] = Z_min + ( Z_max - Z_min ) * (double) i / (double) ns;

	plspal0( "cmap0_black_on_white.pal" );
	plspal1( "cmap1_blue_red.pal", 1 );

	plcol0( 1 );
	plschr( 0., 0.5 );
    plenv( E_min, E_max, N_min, N_max, 1, 0 );

    // Plot the data that was prepared above.
	 for ( i = 1; i <= k; i++ )
    {
		plcol1( ( singleElevation[i] - Z_min ) / ( Z_max - Z_min ) );
		plstring( 1, &singleEasting[i], &singleNorthing[i],"#(729)"  );
	 }
	plschr( 0., 0.5 );

	// draw colorbar
	if ( colorbar1 )
	{
		// Smaller text
		plschr( 0.0, 0.75 );
		// Small ticks on the vertical axis
		plsmaj( 0.0, 0.5 );
		plsmin( 0.0, 0.5 );

		//printf("\nGenerating colorbar");
		num_values1[0] = ns + 1;
		values1[0]     = shedge1;
		plcolorbar( &colorbar_width, &colorbar_height,
			PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL, 0,
			0.005, 0.0, 0.0375, 0.875, 0, 1, 1, 0.0, 0.0,
			cont_color, cont_width,
			n_labels, label_opts, labels,
			n_axis_opts, axis_opts,
			axis_ticks, axis_subticks,
			num_values1, (const double * const *) values1 );

		// Reset text and tick sizes
		plschr( 0.0, 1.0 );
		plsmaj( 0.0, 1.0 );
		plsmin( 0.0, 1.0 );
	}
	plschr( 0.,1.0 );
	plmtex("b",1.7, 0.5, 0.5, "Easting [m]");
	plmtex("l",3.5, 0.5, 0.5, "Northing [m]");
	plmtex("t",1.5, 0.5, 0.5, "LAST: Plan view ");
	save_plot( f_name_last1 );

	////////////////////////////////////////////////////////////////////////
	fclose(out_sum);
	fclose(in_fp2);
	fclose(out_fp2);
	fclose(in_len);

	plend();
    exit( 0 );
	
	return(0);

 } //end of main()

 // save plot
void save_plot( char *fname )
{
    PLINT cur_strm, new_strm;

    plgstrm( &cur_strm );    // get current stream
    plmkstrm( &new_strm );   // create a new one

    plsdev( "psc" );         // new device type. Use a known existing driver
    plsfnam( fname );        // file name

    plcpstrm( cur_strm, 0 ); // copy old stream parameters to new stream
    plreplot();              // do the save
    plend1();                // close new device

    plsstrm( cur_strm );     // and return to previous one
}
	

 

 