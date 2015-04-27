//***********************************************************************
//**               Module 4 of Canopy Height Program                   **
//***********************************************************************
//** Computes Leaf Area Index and Canopy Height Profiles               **
//** in gridded datasets                                               **
//**                                                                   **
//** Authors: Karolina Fieber, Ian Davenport                           **
//** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **
//** University of Reading, UK              November 2014 - March 2015 **
//***********************************************************************

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h> 

/************************************************************************************************************************************************/
/*																FUNCTIONS																		*/
/************************************************************************************************************************************************/
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

int readinstring(FILE *fp, char *string) // to cope with spaces in the file location
{
     char instring[999];
     int i;
     instring[0]=fgetc(fp);
     for (i=1;(instring[i-1]!='\n')&&(i<999);i++)
           instring[i]=fgetc(fp);
     if (i<999)
     {
           instring[i-1]='\0';
           sprintf(string,"%s",instring);
           return (0);
     }
     else
           return(-1);
}

/************************************************************************************************************************************/
/*													P R O G R A M																	*/
/************************************************************************************************************************************/

main()
{
	/*****************/
	/* File pointers */
	/*****************/

	FILE	*in_fp,*in_fp2,*in_fp3,*in_fp4;
	FILE	*out_fp1,*out_fp2,*out_fp3,*out_fp4,*out_fp5,*out_sum,*out_fp8,*out_fp9;
	double	**LAI_array,**refl_array,**mean_single_ground_final,**CHP_array, **LAI_array_corrected;


	
	/*************/
	/* Variables */
	/*************/

	char	in_filename[200], in_filename2[200],in_filename3[200],in_filename4[200];
    char	out_filename1[200],out_filename2[200],out_filename3[200],out_filename4[200],out_filename5[200],out_filename7[200],out_filename8[200],out_filename9[200];
	char	identifier_outX[25], identifier_outY[25],input_YN[2];
	char	characX[50],characY[50], filename[200], *DTM_name;
	int	 	DTM_choice;

	double	plot_reflectance_ratio,Total_LAI_corrected;

	int		length_data,length_dtm;
	int		ygrid_total,xgrid_total;
	int		ygrid,xgrid, count_no_ground, count_no_data;
	double	x_ll,y_ll,ycell_size, xcell_size, ydim, xdim;

	int		count_waveform_grid, count_ground_single_grid; 
	double	single_ground, count_ratio_reflectance;


	int		index, train_counter, number_of_bins, number_of_pulses_param, new_n, pulse_counter;

	double	x[60],x_revers[60],xdata[500];

	int		i,j;
	double	Start_E,Start_N,Start_Z,Vector_E,Vector_N,Vector_Z,E[20],N[20],Z[20];

	double	Z_flat[10],ground_elev, noise;
	int		count,i_ground,begining,ending,begining_bin,offset_bin,last_canopy,last_ground;

	double	ground_E,ground_N,ground_Z,ground_no,ground_empty, mean_ground,begining_Z,canopy_end,ground_end;
	double	energy_veg,energy_ground,mean_single_ground;
	double	canopy_return,ground_return;
	int		zero,last_canopy_bin,last_ground_bin,last_canopy_bin_default;

	double	total_lai,overall_lai,energy_veg_total,energy_ground_total,ratio_refl_total;
	double	Total_LAI,Total_refl,LAI_maximum,chp_total_sum_sum;
	
	double	*amplitude,*amplitude_revers,*amp_noisefree,*area;
	double	*energy,*energy_sum,*height_bin;
	double	*energy_ground_array,*energy_veg_array,*ratio_refl;
	double	*ground_E_array,*ground_N_array,*ground_Z_array,*ground_no_array,*ground_empty_array,*NEW_N1;
	double	*cum_area,*closure,*lai,*chp,*chp_grid;
	double	*lai_plot,*closure_plot,*chp_plot,*chp_plot_grid,*chp_total_sum,*chp_total,*chp_total_norm;
	
	double	total_area,ratio,total, reminder, reminderY;
	double	profile_height, tree_height;
	int		PulseID_param, pulseID,Sampling_unit,Pulse_descriptor,transmitted_samples,TotalNoOfSegments, SegmentNo, segment;
	double	gps_time,waveform_start_easting,waveform_start_northing,waveform_start_elevation;
	double	start_pulse_waveform[100],peak_elev;

	int		sample, profile_bins,first_sample_offset,length_dtm_last;
	int		trans_width; 
	double	width_add,vectorZ;
	double	b, single_amplitude, single_width,amplitude_trans;

	DTM_name=(char *)calloc(15,sizeof(char));

	x_ll=y_ll=xcell_size=ycell_size=plot_reflectance_ratio=-999;
	sample=trans_width=-999;

	printf("***********************************************************************\n");
	printf("**               Module 4 of Canopy Height Program                   **\n");
	printf("***********************************************************************\n");
	printf("** Computes Leaf Area Index and Canopy Height Profiles               **\n");
	printf("** in gridded datasets                                               **\n");
	printf("**                                                                   **\n");
	printf("** Authors: Karolina Fieber, Ian Davenport                           **\n");
	printf("** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	printf("** University of Reading, UK                              March 2015 **\n");
	printf("***********************************************************************\n\n");

	///////////////////
	// USER INPUT
	///////////////////
	printf("Enter the name of the input file (without '_fulwvs' suffix): ");
	scanf("%s",&filename);

	//OPENING FILES for writing
	sprintf(out_filename7,"%s_info_Module4.txt", filename);
	out_sum=fopen(out_filename7,"w");
	if (!out_sum)
	{
		printf("Can't open output data file <%s>\n",out_filename7);
 		return(-1);
	}

	fprintf(out_sum,"***********************************************************************\n");
	fprintf(out_sum,"**               Module 4 of Canopy Height Program                   **\n");
	fprintf(out_sum,"***********************************************************************\n");
	fprintf(out_sum,"** Computes Leaf Area Index and Canopy Height Profiles               **\n");
	fprintf(out_sum,"** in gridded datasets                                               **\n");
	fprintf(out_sum,"**                                                                   **\n");
	fprintf(out_sum,"** Authors: Karolina Fieber, Ian Davenport                           **\n");
	fprintf(out_sum,"** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	fprintf(out_sum,"** University of Reading, UK                              March 2015 **\n");
	fprintf(out_sum,"***********************************************************************\n\n");

	/* Reading from files*/
	sprintf(in_filename,"%s_fulwvs.txt",filename);
	in_fp=fopen(in_filename,"r");
	if (!in_fp)
	{
		printf("\nCan't open input data file <%s>\n",in_filename);
		fprintf(out_sum,"\nCan't open input data file <%s>\n",in_filename);
 		return(-1);
	}
	else
	{
		printf("\n\tReading from file <<%s>>\n",in_filename);
		fprintf(out_sum,"\n\tReading from file <<%s>>\n",in_filename);
	}

	sprintf(in_filename2,"%s_received_param.txt", filename);
	in_fp2=fopen(in_filename2,"r");
	if (!in_fp2)
	{
		printf("Can't open input data file <%s>\n",in_filename2);
		fprintf(out_sum,"Can't open input data file <%s>\n",in_filename);
 		return(-1);
	}
	else
	{
		printf("\tReading from file <<%s>>\n",in_filename2);
		fprintf(out_sum,"\tReading from file <<%s>>\n",in_filename2);
	}
	
	sprintf(in_filename3,"%s_param.txt",filename);
	in_fp3=fopen(in_filename3,"r");
	if (!in_fp3)
	{
		printf("Can't open input data file <%s>\n",in_filename3);
		fprintf(out_sum,"Can't open input data file <%s>\n",in_filename);
 		return(-1);
	}
	else
	{
		printf("\tReading from file <<%s>>\n",in_filename3);
		fprintf(out_sum,"\tReading from file <<%s>>\n",in_filename3);
	}

	// Read in the length file
	length_data=readininteger(in_fp3); // total number of waveforms
	length_dtm=readininteger(in_fp3); // number of single returns
	length_dtm_last=readininteger(in_fp3); //total number of points	detected
	sample=readininteger(in_fp3); // sampling unit

	amplitude_trans=readindouble(in_fp3); // amplitude of transmitted pulse
	width_add=readindouble(in_fp3); //  width of transmitted pulse
	if (width_add-(int)width_add>0.51)
		trans_width=(int)width_add+1;
	else
		trans_width=(int)width_add;

	single_amplitude=readindouble(in_fp3); // mean amplitude of single returns
	single_width=readindouble(in_fp3); // mean width of single returns

	noise=readindouble(in_fp3); // mean noise in transmitted waveform
	vectorZ=readindouble(in_fp3); // elevation vector

	b=(vectorZ*100-(int)(vectorZ*100));

	if ((fabs(vectorZ*100-(int)(vectorZ*100)))>0.51)
		vectorZ=(fabs((int)(vectorZ*100))+1)/100;
	else
		vectorZ=(fabs((int)(vectorZ*100))+0.5)/100;

	DTM_choice=1;
	printf("\nWhich file with ground returns would you like to select \n\tsingle return = 1 (default); \n\tlast return = 2;\n\tyour own = 3? ");
	scanf("%d",&DTM_choice);
	while (DTM_choice!=1 && DTM_choice!=2 && DTM_choice!=3 )
	{
		printf("\n\tWrong input: ");
		scanf("%d",&DTM_choice);
	}

	if (DTM_choice==1)
	{
		DTM_name="d_single.txt";
	}
	else if (DTM_choice==2)
	{
		length_dtm=length_data;
		DTM_name="d_last.txt";;
	}
	else if (DTM_choice==3)
	{
		length_dtm=length_dtm_last;
		DTM_name="d_own.txt";;
	}


	sprintf(in_filename4,"%s_groun%s",filename, DTM_name);
	in_fp4=fopen(in_filename4,"r");
	if (!in_fp4)
	{
		printf("Can't open input data file <%s>\n",in_filename4);
		fprintf(out_sum,"Can't open input data file <%s>\n",in_filename);
 		return(-1);
	}
	else
	{
		printf("\n\tReading from file <<%s>>\n",in_filename4);
		fprintf(out_sum,"\tReading from file <<%s>>\n",in_filename4);
	}

	fprintf(out_sum,"\n***********************************************************************\n");

	/// Input parameters
	printf("\nEnter lower left corner of the area of interest (Easting) [m]: ");
	scanf("%lf",&x_ll);
	printf("Enter lower left corner of the area of interest (Northing) [m]: ");
	scanf("%lf",&y_ll);

	printf("\n\nEnter the grid cell size X (Easting/columns) [m]: ");
	scanf("%lf",&xcell_size);
	printf("Enter the grid cell size Y (Northing/rows) [m]: ");
	scanf("%lf",&ycell_size);

	
	printf("\n\nEnter site dimension X (Easting/columns) [m]: ");
	scanf("%lf",&xdim);
	reminder=xdim/xcell_size;
	while (reminder-(int)reminder!=0)
	{
		printf("\n\nEnter a different site dimension X (Easting/columns) that would be a multiplication of grid cell size X [m]: ");
		scanf("%lf",&xdim);
		reminder=xdim/xcell_size;
	}
	
	printf("Enter site dimension Y (Northing/rows) [m]: ");
	scanf("%lf",&ydim);
	reminderY=ydim/ycell_size;
	while (reminderY-(int)reminderY!=0)
	{
		printf("\n\nEnter a different site dimension Y (Northing/rows) that would be a multiplication of grid cell size Y [m]: ");
		scanf("%lf",&ydim);
		reminderY=ydim/ycell_size;
	}
	
	printf("\nEnter approximate maximum height of the trees [m]: ");
	scanf("%lf",&tree_height);

	
	printf("\nEnter default reflectance ratio \n\t0.5 for 1550nm wavelength, \n\t2.0 for 1064nm wavelength: ");
	scanf("%lf",&plot_reflectance_ratio);

	
	//to output the gridding size to the name of the output files
	sprintf(characX,"%0.0lf",xcell_size);
	sprintf(characY,"%0.0lf",ycell_size);
	
	ygrid_total=(int)(ydim/ycell_size); /*number of grid cell -rows*/
	xgrid_total=(int)(xdim/xcell_size); /*number of grid cell -columns*/
	sprintf(identifier_outX,"%s", characX);
	sprintf(identifier_outY,"%s", characY);

	printf("\nSpecify separation of ground/vegetation in number of bins \n\t(if you want it to be calculated by the program: -999): ");
	scanf("%d",&last_canopy_bin_default);
	while ((last_canopy_bin_default>=50 || last_canopy_bin_default<=0) && last_canopy_bin_default!=-999)
	{
		printf("\n\nWrong input");
		scanf("%lf",&last_canopy_bin_default);
	}

	profile_bins=(int)(tree_height/(sample*vectorZ))+(int)(14/sample);// adding extra about 2m in case the trees are taller than specified
	profile_height=profile_bins*sample*vectorZ;

	printf("\nWould you like to generate files of Energy, Closure etc. profiles in grid? ");
	scanf("%s",&input_YN);
	while ((input_YN[0]!='Y'&&input_YN[0]!='y')&&(input_YN[0]!='N'&&input_YN[0]!='n'))
	{
		printf("\nWRONG INPUT: Y/N? ");
		scanf("%s",&input_YN);
	}

	mean_ground=-999;

	// END OF USER INPUT

	energy_ground_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	energy_veg_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	
	ratio_refl=(double *)calloc(length_dtm_last+2,sizeof(double));
	ground_E_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	ground_N_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	ground_Z_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	ground_no_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	ground_empty_array=(double *)calloc(length_dtm_last+2,sizeof(double));
	NEW_N1=(double *)calloc(length_dtm_last+2,sizeof(double));

	amplitude=(double *)calloc(5000,sizeof(double));
	amplitude_revers=(double *)calloc(5000,sizeof(double));
	amp_noisefree=(double *)calloc(5000,sizeof(double));
	area=(double *)calloc(5000,sizeof(double));

	energy=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	energy_sum=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	height_bin=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	cum_area=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	closure=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	closure_plot=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	lai=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	lai_plot=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp_plot=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp_grid=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp_plot_grid=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp_total_sum=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp_total=(double *)calloc(profile_bins+3*trans_width,sizeof(double));
	chp_total_norm=(double *)calloc(profile_bins+3*trans_width,sizeof(double));

	mean_single_ground_final=(double **)calloc(ygrid_total+2,sizeof(double));
	for (j=0;j<=ygrid_total+1;j++)
		mean_single_ground_final[j]=(double *)calloc(xgrid_total+1,sizeof(double));

	for (i=1;i<=ygrid_total;i++)
		for (j=1;j<=xgrid_total;j++)
		  {
					mean_single_ground_final[i][j]=(double)(0);
		  }


	LAI_array=(double **)calloc(ygrid_total+2,sizeof(double));
	for (j=0;j<=ygrid_total+1;j++)
		LAI_array[j]=(double *)calloc(xgrid_total+1,sizeof(double));

	for (i=1;i<=ygrid_total;i++)
		for (j=1;j<=xgrid_total;j++)
		  {
					LAI_array[i][j]=(double)(0);
		  }

	LAI_array_corrected=(double **)calloc(ygrid_total+2,sizeof(double));
	for (j=0;j<=ygrid_total+1;j++)
		LAI_array_corrected[j]=(double *)calloc(xgrid_total+1,sizeof(double));

	for (i=1;i<=ygrid_total;i++)
		for (j=1;j<=xgrid_total;j++)
		  {
					LAI_array_corrected[i][j]=(double)(0);
		  }

	refl_array=(double **)calloc(ygrid_total+2,sizeof(double));
	for (j=0;j<=ygrid_total+1;j++)
		refl_array[j]=(double *)calloc(xgrid_total+1,sizeof(double));

	for (i=1;i<=ygrid_total;i++)
		for (j=1;j<=xgrid_total;j++)
		  {
					refl_array[i][j]=(double)(0);
		  }

	CHP_array=(double **)calloc(ygrid_total*xgrid_total+2,sizeof(double));
	for (j=0;j<=ygrid_total*xgrid_total+1;j++)
	{
		CHP_array[j]=(double *)calloc(profile_bins+3*trans_width+1,sizeof(double));
	}

	for (i=1;i<=ygrid_total*xgrid_total;i++)
	{
		for (j=1;j<=profile_bins+3*trans_width;j++)
		  {
					CHP_array[i][j]=(double)(0);
		  }
	}
	
	printf("\n*************************************************************************\n");
	fprintf(out_sum,"\nUSER INPUT VALUES\n");
	fprintf(out_sum,"\nThe lower left corner coordinates are: \t%.2lfm, %.2lfm", x_ll,y_ll);
	fprintf(out_sum,"\nThe grid cell size is: \t\t\t%.2lfm by %.2lfm", xcell_size,ycell_size);
	fprintf(out_sum,"\nThe site dimensions are: \t\t%.2lfm by %.2lfm (%.2lfm2) \n\tmaking total number of cells \t%.0lf (%.0lf by %.0lf)\n\n", xdim,ydim,xdim*ydim, ((xdim/xcell_size)*(ydim/ycell_size)),xdim/xcell_size,ydim/ycell_size);
	fprintf(out_sum,"Sampling vertical distance is: \t%dns \n\twhich equals to approx: %.3fm\n\n", sample,sample*vectorZ);
	fprintf(out_sum, "Transmitted pulse width is: \t%d\n", trans_width);
	fprintf(out_sum, "Default reflectance ratio is: \t%0.2lf\n", plot_reflectance_ratio);
	fprintf(out_sum,"\nTree height specified by user is \t%0.2lfm\n", tree_height);
	fprintf(out_sum,"Profile height length is: \t\t%.2lfm \n\teqating to %d bins of %.2lfcm\n", profile_height, profile_bins,sample*vectorZ);
	if (last_canopy_bin_default!=-999)
		fprintf(out_sum,"\nDefault location where the canopy ends: %dth bin AGL, \n\tequating to %0.2lfm AGL\n", last_canopy_bin_default,last_canopy_bin_default*sample*vectorZ);
	else
		fprintf(out_sum,"\nNo default location where the canopy ends specified by user (-999).\n");
	fprintf(out_sum,"\nWould you like to generate files of Energy, Closure etc. profiles in grid? %s\n", input_YN);

	/*files OUT*/

	if (input_YN[0]=='Y'||input_YN[0]=='y')
	{
		sprintf(out_filename1,"%s_%smx%sm_1_Energy_grid.txt",filename,identifier_outX, identifier_outY);
		out_fp1=fopen(out_filename1,"w");
		if (!out_fp1)
		{
			printf("Can't open output data file <%s>\n",out_filename1);
 			return(-1);
		}
		fprintf(out_sum,"\n\tWriting to file <<%s>>", out_filename1);
		printf("\n\tWriting to file <<%s>>", out_filename1);


		sprintf(out_filename2,"%s_%smx%sm_2_Closure_grid.txt", filename,identifier_outX, identifier_outY);
		out_fp2=fopen(out_filename2,"w");
		if (!out_fp2)
		{
			printf("Can't open output data file <%s>\n",out_filename2);
 			return(-1);
		}
		fprintf(out_sum,"\n\tWriting to file <<%s>>", out_filename2);
		printf("\n\tWriting to file <<%s>>", out_filename2);

		sprintf(out_filename3,"%s_%smx%sm_3_LAI_grid.txt", filename,identifier_outX, identifier_outY);
		out_fp3=fopen(out_filename3,"w");
		if (!out_fp3)
		{
			printf("Can't open output data file <%s>\n",out_filename3);
 			return(-1);
		}
		fprintf(out_sum,"\n\tWriting to file <<%s>>", out_filename3);
		printf("\n\tWriting to file <<%s>>", out_filename3);

		sprintf(out_filename4,"%s_%smx%sm_4_CHP_grid.txt", filename,identifier_outX, identifier_outY);
		out_fp4=fopen(out_filename4,"w");
		if (!out_fp4)
		{
			printf("Can't open output data file <%s>\n",out_filename4);
 			return(-1);
		}
		fprintf(out_sum,"\n\tWriting to file <<%s>>", out_filename4);
		printf("\n\tWriting to file <<%s>>", out_filename4);

		sprintf(out_filename5,"%s_%smx%sm_5_Reflectance.txt", filename,identifier_outX, identifier_outY);
		out_fp5=fopen(out_filename5,"w");
		if (!out_fp5)
		{
			printf("Can't open output data file <%s>\n",out_filename5);
 			return(-1);
		}
		fprintf(out_sum,"\n\tWriting to file <<%s>>", out_filename5);
		printf("\n\tWriting to file <<%s>>", out_filename5);
	}

	sprintf(out_filename8,"%s_%smx%sm_LAI_final.txt", filename,identifier_outX, identifier_outY);
	out_fp8=fopen(out_filename8,"w");
	if (!out_fp8)
	{
		printf("Can't open output data file <%s>\n",out_filename8);
 		return(-1);
	}
	fprintf(out_sum,"\n\tWriting to file <<%s>>", out_filename8);
	printf("\n\tWriting to file <<%s>>", out_filename8);

	sprintf(out_filename9,"%s_%smx%sm_CHP_final.txt", filename,identifier_outX, identifier_outY);
	out_fp9=fopen(out_filename9,"w");
	if (!out_fp9)
	{
		printf("Can't open output data file <%s>\n",out_filename9);
 		return(-1);
	}
	fprintf(out_sum,"\n\tWriting to file <<%s>>\n", out_filename9);
	printf("\n\tWriting to file <<%s>>\n", out_filename9);
	
	fprintf(out_sum,"\n*************************************************************************\n");
	fprintf(out_sum,"\nPROGRAM OUTPUT VALUES\n");
	
	for (i=1;i<=profile_bins+3*trans_width;i++)
	{
		height_bin[i]=profile_height-(i-1)*sample*vectorZ;
	}

	if (input_YN[0]=='Y'||input_YN[0]=='y')
	{
		fprintf(out_fp1,"Cell no,,,");
		fprintf(out_fp2,"Cell no,,");
		fprintf(out_fp3,"Cell no,,,");
		fprintf(out_fp4,"Cell no,,,");
		fprintf(out_fp5,"Cell no: Vegetation Energy, Ground Energy, Reflectance Ratio, Mean single ground energy\n");
	}
	//fprintf(out_fp9,"Cell no,");

	fprintf(out_fp9,"%d\n",profile_bins+3*trans_width-1);

	for (i=1;i<=profile_bins+3*trans_width-1;i++)
	{
		if (input_YN[0]=='Y'||input_YN[0]=='y')
		{
			fprintf(out_fp1,"%f,",height_bin[i]-(sample*vectorZ/2));
			fprintf(out_fp2,"%f,",height_bin[i]-(sample*vectorZ/2));
			fprintf(out_fp3,"%f,",height_bin[i]-(sample*vectorZ/2));
			fprintf(out_fp4,"%f,",height_bin[i]-(sample*vectorZ/2));
		}
		fprintf(out_fp9,"%f,",height_bin[i]-(sample*vectorZ/2));
	}
	if (input_YN[0]=='Y'||input_YN[0]=='y')
	{
		fprintf(out_fp1,"\n");
		fprintf(out_fp2,"\n");
		fprintf(out_fp3,"\n");
		fprintf(out_fp4,"\n");
	}
	fprintf(out_fp9,"\n");
//END OF OPENING FILES AND CREATING OUTPUT FILES	

	/*****************/
	/* P R O G R A M */
	/*****************/

	pulse_counter=0;
	count_no_ground=0;// number of cells with no ground return
	count_no_data=0; // number of cell with no data in

	{
		for (i_ground=0;i_ground<=length_dtm;i_ground++)
		{
			ground_no_array[i_ground]=0;
			ground_E_array[i_ground]=0;
			ground_N_array[i_ground]=0;
			ground_Z_array[i_ground]=0;
			ground_empty_array[i_ground]=0;
		}

		for (i_ground=1;i_ground<=length_dtm;i_ground++)
		{
			ground_no_array[i_ground]=readininteger(in_fp4);
			ground_E_array[i_ground]=readindouble(in_fp4);
			ground_N_array[i_ground]=readindouble(in_fp4);
			ground_Z_array[i_ground]=readindouble(in_fp4);
			ground_empty_array[i_ground]=readininteger(in_fp4);
		}
		
		//For each grid cell
		for (ygrid=0;ygrid<=ygrid_total-1;ygrid++) /*grid definition - row*/
		{
			for (xgrid=0;xgrid<=xgrid_total-1;xgrid++) /* column*/
			{
				count_waveform_grid=0; /*count how many waveforms are in each grid square*/
				count_ground_single_grid=0; /*count how single ground returns are in each grid square*/
				single_ground=0;
				mean_single_ground=0; 
				count_ratio_reflectance=0;
				for (i=0; i<=profile_bins+3*trans_width-1; i++)
					NEW_N1[i]=0;

				rewind(in_fp3);
				rewind(in_fp2);
				rewind(in_fp);

				printf("\n Cell: [%d,%d]",ygrid+1, xgrid+1);
				

				for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
					energy_sum[j]=0;
				
				
				for (index=1;index<=length_data;index++) 
				{
					pulse_counter++;
					/*printf("%d,", pulse_counter);*/

					//if (pulse_counter==19191)
						//printf("hi");


					for (j=0;j<=5000-1;j++) /*zero area array*/
						area[j]=0;
					for (j=0;j<=5000-1;j++) /*this is to make sure that long waveform trains of previous waveform is erased.*/
						amplitude[j]=0;
					for (j=0;j<=5000-1;j++) /*this is to make sure that long waveform trains of previous waveform is erased.*/
						amplitude_revers[j]=0;

					/*{For all waveforms:... 
					... 1. Read in data 
					... 2. Check which ones fall in particular block
					...     a. calculate mean ground level
					...     b. Subtract mean noise from each waveform 
					...     c. Find the begining of the waveform
					...     d. Find the end of the last pulse
					...     e. Calculate the area of each bin and align them 
					...     f. Calulate the mean energy of single peak ground returns */
					
					/**************************************/
					/*1.Read in data*/
					/* read amplitude train*/

					pulseID=readininteger(in_fp); //pulseID
					gps_time=readindouble(in_fp); //GPS time

					Sampling_unit=readininteger(in_fp);					// Sampling unit eg. 1ns
					Pulse_descriptor=readininteger(in_fp);				// Pulse descriptor in the original file
					Start_E=readindouble(in_fp);			// anchor location easting
					Start_N=readindouble(in_fp);			// anchor location northing
					Start_Z=readindouble(in_fp);			// anchor location elevation
					transmitted_samples=readininteger(in_fp);			// number of samples in transmitted waveform
					for (train_counter=0;train_counter<=transmitted_samples-1;train_counter++)
					{
						start_pulse_waveform[train_counter]=readindouble(in_fp);
					}
				
					Vector_E=readindouble(in_fp);		// vector easting
					Vector_N=readindouble(in_fp);		// vector northing
					Vector_Z=readindouble(in_fp);		// vector elevation
					TotalNoOfSegments=readininteger(in_fp);				// number of segments of received waveform

					for (segment=1;segment<=TotalNoOfSegments;segment++)
					{
						SegmentNo=readininteger(in_fp);				//segment number
						first_sample_offset=readininteger(in_fp);		//offset to the first sample of that segment

						waveform_start_easting=readindouble(in_fp);	//easting coordinate of the first sample in the segment
						waveform_start_northing=readindouble(in_fp);	//northing coordinate of the first sample in the segment
						waveform_start_elevation=readindouble(in_fp);	//elevation coordinate of the first sample in the segment
						number_of_bins=readininteger(in_fp);					//number of samples in received waveform (particular segment)
						for (train_counter=0;train_counter<=number_of_bins-1;train_counter++)
						{
							amplitude[train_counter]=readindouble(in_fp);
						}
				
						peak_elev=readindouble(in_fp);					// elevation of the bin with the highest amplitude
						
						//read in peak parameters
						xdata[0]=0;
						PulseID_param=readininteger(in_fp2);

						//if (PulseID_param!=pulseID)
							//printf("Trouble, %d vs %d\n",pulseID,PulseID_param);

						if (PulseID_param!=-999) // if the waveform was disregarded or not
						{
							number_of_pulses_param=readininteger(in_fp2);	
							for (train_counter=1;train_counter<=number_of_bins-1; train_counter++)
							{
								xdata[train_counter]=train_counter-1;
							}

							/*read peak parameters*/
					
							for (j=0;j<=59;j++) /*to zero the array*/
								x[j]=0;
							
							for (j=0;j<=59;j++) /*to zero the array*/
								x_revers[j]=0;

							for (i=1;i<=number_of_pulses_param;i++)   
							{
									x[(i-1)*3 + 1] = readindouble(in_fp2); /*peak position x0(1),x0(4),x0(7)...*/
									x[(i-1)*3 + 2] = readindouble(in_fp2); /*amplitude x0(2),x0(5),x0(8)...*/
									x[(i-1)*3 + 3] = readindouble(in_fp2); /*standard deviation x0(3),x0(6),x0(9)...*/
							}

							// check for 'upside down' recording
							if (Vector_Z>0 &&Start_Z<=waveform_start_elevation)
							{
								//new start pulse position as the 'last' sample of the 'upside down' waveform
								Start_E=Start_E+Vector_E*(first_sample_offset+number_of_bins);
								Start_N=Start_N+Vector_N*(first_sample_offset+number_of_bins);
								Start_Z=Start_Z+Vector_Z*(first_sample_offset+number_of_bins);
								first_sample_offset=0;
								for (train_counter=0;train_counter<=number_of_bins-1;train_counter++)
								{
									amplitude_revers[number_of_bins-1-train_counter]=amplitude[train_counter];
								}
								
								for (train_counter=0;train_counter<=number_of_bins-1;train_counter++)
								{
									amplitude[train_counter]=amplitude_revers[train_counter];
								}

								for (i=1;i<=number_of_pulses_param;i++)   
								{
									x_revers[(number_of_pulses_param-i)*3 + 1] =number_of_bins-x[(i-1)*3 + 1];
									x_revers[(number_of_pulses_param-i)*3 + 2] =x[(i-1)*3 + 2];
									x_revers[(number_of_pulses_param-i)*3 + 3] =x[(i-1)*3 + 3];
								}

								for (i=1;i<=number_of_pulses_param;i++)   
								{
									x[(i-1)*3 + 1]=x_revers[(i-1)*3 + 1];
									x[(i-1)*3 + 2]=x_revers[(i-1)*3 + 2];
									x[(i-1)*3 + 3]=x_revers[(i-1)*3 + 3];
								}
								Vector_Z=-Vector_Z;
								Vector_N=-Vector_N;
								Vector_E=-Vector_E;
							}
				
							for (j=0;j<=19;j++) /*zero*/
							{
								E[j]=0;
								N[j]=0;
								Z[j]=0;
							}

							for (i=1;i<=number_of_pulses_param;i++) 
							{
								E[i]=Start_E+Vector_E*(first_sample_offset+x[(i-1)*3+1]);   /*easting*/
								N[i]=Start_N+Vector_N*(first_sample_offset+x[(i-1)*3+1]);   /*northing*/
								Z[i]=Start_Z+Vector_Z*(first_sample_offset+x[(i-1)*3+1]);   /*elevation*/
							}
							new_n=number_of_pulses_param;

						/************* end of reading*******************/

						/*2. check which waveforms fall in particular block*/
							if ((E[number_of_pulses_param]>=(x_ll+xgrid*xcell_size)&&E[number_of_pulses_param]<=(x_ll+(xgrid+1)*xcell_size))&&(N[number_of_pulses_param]>=(y_ll+ygrid*ycell_size)&&N[number_of_pulses_param]<=y_ll+(ygrid+1)*ycell_size))
							{
								for (j=0;j<=9;j++) 
										Z_flat[j]=0;

								ground_elev=0;
								count=0;
                
								/*2a. calculate mean local ground level for the last pulse within waveform*/
								for (i_ground=1;i_ground<=length_dtm;i_ground++)
								{
									ground_no=ground_no_array[i_ground];
									ground_E=ground_E_array[i_ground];
									ground_N=ground_N_array[i_ground];
									ground_Z=ground_Z_array[i_ground];
									ground_empty=ground_empty_array[i_ground];
						
									if (((ground_E<=E[number_of_pulses_param]+0.5)&&ground_E>=(E[number_of_pulses_param]-0.5))&&(ground_N<=(N[number_of_pulses_param]+0.5)&&ground_N>=(N[number_of_pulses_param]-0.5)))
									{
										ground_elev=ground_elev+ground_Z;
										count=count+1;
										mean_ground=ground_elev/count;
									}
								}

								if (count<1) /*if the there are less than one points of dtm within the proximity of the waveform check larger area*/
								{
									for (i_ground=1;i_ground<=length_dtm;i_ground++)
									{
										ground_no=ground_no_array[i_ground];
										ground_E=ground_E_array[i_ground];
										ground_N=ground_N_array[i_ground];
										ground_Z=ground_Z_array[i_ground];
										ground_empty=ground_empty_array[i_ground];
						
										if (((ground_E<=E[number_of_pulses_param]+1.0)&&ground_E>=(E[number_of_pulses_param]-1.0))&&(ground_N<=(N[number_of_pulses_param]+1.0)&&ground_N>=(N[number_of_pulses_param]-1.0)))
										{
											ground_elev=ground_elev+ground_Z;
											count=count+1;
											mean_ground=ground_elev/count;
										}
									}
								}

								if (count<1) /*if the there are less than one points of dtm within the proximity of the waveform check larger area*/
								{
									for (i_ground=1;i_ground<=length_dtm;i_ground++)
									{
										ground_no=ground_no_array[i_ground];
										ground_E=ground_E_array[i_ground];
										ground_N=ground_N_array[i_ground];
										ground_Z=ground_Z_array[i_ground];
										ground_empty=ground_empty_array[i_ground];
						
										if (((ground_E<=E[number_of_pulses_param]+1.5)&&ground_E>=(E[number_of_pulses_param]-1.5))&&(ground_N<=(N[number_of_pulses_param]+1.5)&&ground_N>=(N[number_of_pulses_param]-1.5)))
										{
											ground_elev=ground_elev+ground_Z;
											count=count+1;
											mean_ground=ground_elev/count;
										}
									}
								}

								if (count<1) /*if the there are less than one points of dtm within the proximity of the waveform check larger area*/
								{
									for (i_ground=1;i_ground<=length_dtm;i_ground++)
									{
										ground_no=ground_no_array[i_ground];
										ground_E=ground_E_array[i_ground];
										ground_N=ground_N_array[i_ground];
										ground_Z=ground_Z_array[i_ground];
										ground_empty=ground_empty_array[i_ground];
						
										if (((ground_E<=E[number_of_pulses_param]+2.5)&&ground_E>=(E[number_of_pulses_param]-2.5))&&(ground_N<=(N[number_of_pulses_param]+2.5)&&ground_N>=(N[number_of_pulses_param]-2.5)))
										{
											ground_elev=ground_elev+ground_Z;
											count=count+1;
											mean_ground=ground_elev/count;
										}
									}
								}

								if (count<1) /*if the there are less than one points of dtm within the proximity of the waveform check larger area*/
								{
									for (i_ground=1;i_ground<=length_dtm;i_ground++)
									{
										ground_no=ground_no_array[i_ground];
										ground_E=ground_E_array[i_ground];
										ground_N=ground_N_array[i_ground];
										ground_Z=ground_Z_array[i_ground];
										ground_empty=ground_empty_array[i_ground];
						
										if (((ground_E<=E[number_of_pulses_param]+3.5)&&ground_E>=(E[number_of_pulses_param]-3.5))&&(ground_N<=(N[number_of_pulses_param]+3.5)&&ground_N>=(N[number_of_pulses_param]-3.5)))
										{
											ground_elev=ground_elev+ground_Z;
											count=count+1;
											mean_ground=ground_elev/count;
										}
									}
								}

								new_n=number_of_pulses_param;
								Z_flat[new_n]=Z[number_of_pulses_param]-mean_ground;

								/*to check if there is something way below ground level - upto two
								% ringing echoes; if it is new mean ground level is
								% calculated */
								if (Z_flat[new_n]<-0.5&&new_n>1)
								{	 
									new_n=number_of_pulses_param-1;
										E[new_n]=Start_E+Vector_E*(first_sample_offset+x[(number_of_pulses_param-2)*3+1]);   /*easting*/
										N[new_n]=Start_N+Vector_N*(first_sample_offset+x[(number_of_pulses_param-2)*3+1]);   /*northing*/
										Z[new_n]=Start_Z+Vector_Z*(first_sample_offset+x[(number_of_pulses_param-2)*3+1]);   /*elevation*/

										for (j=0;j<=9;j++) /*zero*/
										Z_flat[j]=0;

										ground_elev=0;
										count=0;

										for (i_ground=1;i_ground<=length_dtm;i_ground++) /*new mean ground level calculation*/
										{
										ground_no=ground_no_array[i_ground];
										ground_E=ground_E_array[i_ground];
										ground_N=ground_N_array[i_ground];
										ground_Z=ground_Z_array[i_ground];
										ground_empty=ground_empty_array[i_ground];

										if (((ground_E<=E[number_of_pulses_param]+2)&&ground_E>=(E[number_of_pulses_param]-2))&&(ground_N<=(N[number_of_pulses_param]+2)&&ground_N>=(N[number_of_pulses_param]-2)))
										{
											ground_elev=ground_elev+ground_Z;
											count=count+1;
											mean_ground=ground_elev/count;
										}
										}

									if (count<1)
									{
										for (i_ground=1;i_ground<=length_dtm;i_ground++)
										{
											ground_no=ground_no_array[i_ground];
											ground_E=ground_E_array[i_ground];
											ground_N=ground_N_array[i_ground];
											ground_Z=ground_Z_array[i_ground];
											ground_empty=ground_empty_array[i_ground];
									
											if (((ground_E<=E[number_of_pulses_param]+4)&&ground_E>=(E[number_of_pulses_param]-4))&&(ground_N<=(N[number_of_pulses_param]+4)&&ground_N>=(N[number_of_pulses_param]-4)))
											{	
												ground_elev=ground_elev+ground_Z;
												count=count+1;
												mean_ground=ground_elev/count;
											}
										}
									}
									Z_flat[new_n]=Z[new_n]-mean_ground;
								}

								if (Z_flat[new_n]<-0.5&&new_n>1) /*check for second ringing echoe*/
								{	
										new_n=number_of_pulses_param-1;
										E[new_n]=Start_E+Vector_E*(first_sample_offset+x[(number_of_pulses_param-2)*3+1]);   /*easting*/
										N[new_n]=Start_N+Vector_N*(first_sample_offset+x[(number_of_pulses_param-2)*3+1]);   /*northing*/
										Z[new_n]=Start_Z+Vector_Z*(first_sample_offset+x[(number_of_pulses_param-2)*3+1]);   /*elevation*/

										for (j=0;j<=9;j++) 
										Z_flat[j]=0;

										ground_elev=0;
										count=0;

										for (i_ground=1;i_ground<=length_dtm;i_ground++) /*new mean ground level calculation*/
										{
										ground_no=ground_no_array[i_ground];
										ground_E=ground_E_array[i_ground];
										ground_N=ground_N_array[i_ground];
										ground_Z=ground_Z_array[i_ground];
										ground_empty=ground_empty_array[i_ground];

										if (((ground_E<=E[number_of_pulses_param]+2)&&ground_E>=(E[number_of_pulses_param]-2))&&(ground_N<=(N[number_of_pulses_param]+2)&&ground_N>=(N[number_of_pulses_param]-2)))
										{
											ground_elev=ground_elev+ground_Z;
											count=count+1;
											mean_ground=ground_elev/count;
										}
										}

									if (count<2)
									{
										for (i_ground=1;i_ground<=length_dtm;i_ground++)
										{
											ground_no=ground_no_array[i_ground];
											ground_E=ground_E_array[i_ground];
											ground_N=ground_N_array[i_ground];
											ground_Z=ground_Z_array[i_ground];
											ground_empty=ground_empty_array[i_ground];
									
											if (((ground_E<=E[number_of_pulses_param]+4)&&ground_E>=(E[number_of_pulses_param]-4))&&(ground_N<=(N[number_of_pulses_param]+4)&&ground_N>=(N[number_of_pulses_param]-4)))
											{	
												ground_elev=ground_elev+ground_Z;
												count=count+1;
												mean_ground=ground_elev/count;
											}
										}
									}
									Z_flat[new_n]=Z[new_n]-mean_ground;
								}

								if (mean_ground==-999)
								{
									fprintf(out_sum, "\nWARNING! No ground points available in the area to enable calculation of the ground level!");
									return(-1);
								}
                
								/*end of mean ground level calculation*/
								/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

								/*2b. Subtract mean noise from each waveform*/
								for (j=0;j<=5000-1;j++) 
									amp_noisefree[j]=0;

								for (train_counter=1;train_counter<=number_of_bins; train_counter++)
								{
									amp_noisefree[train_counter]=amplitude[train_counter]-noise;
									if (amp_noisefree[train_counter]<0)
									{
										amp_noisefree[train_counter]=0;
									}
								}

								for (i=(int)(x[(new_n-1)*3+1]+10+1);i<=number_of_bins;i++) /*x value is for amplitude+1 (wartosc 15 binu jest 16ta w kolejce)*/
									amp_noisefree[train_counter]=0;
						
					
								/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
                
								/*2c. Find the begining of the first pulse*/
								if (x[1]>7)
								{
									begining=-999;
									for (i=(int)(x[1]-trans_width*1.5);i<=(int)(x[1]+0.5);i++)
									{
										if ((amp_noisefree[i]>(noise)/2)&&begining<0)
										{
											begining=(int)xdata[i];
										}
									}

									if (begining==-999)
									{
										for (i=(int)(x[1]-trans_width*1.5);i<=(int)(x[1]+1.5);i++)
										{
											if ((amp_noisefree[i]>(noise)/4)&&begining<0)
											{
												begining=(int)xdata[i];
											}
										}
									}

									if (begining==-999 && x[2]<0.002)
									{
										for (i=(int)((x[4])-trans_width*1.5);i<=(int)(x[4]+1.5);i++)
										{
											if ((amp_noisefree[i]>(noise)/4)&&begining<0)
											{
												begining=(int)xdata[i];
											}
										}
									}

									if (begining==-999)
									{
										for (i=1; i<=(int)(x[1]+1.5); i++)
										{
											if ((amp_noisefree[i]>0)&&begining<0)
											{
												begining=(int)xdata[i];
											}
										}
									}
                     
									begining_Z=(Start_Z+Vector_Z*(first_sample_offset+begining))-mean_ground; //relative elevation of the begining of the pulse
								}
                     
								else
								{
									begining=-999;
									for (i=(int)(x[1]-trans_width);i<=(int)(x[1]+0.5);i++)
									{
										if ((amp_noisefree[i]>(noise)/2)&&begining<0)
										{
											begining=(int)xdata[i];
										}
									}
									begining_Z=(Start_Z+Vector_Z*(first_sample_offset+begining))-mean_ground;
								}

								/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
                
								/*2d. Find the end of the last pulse*/
								ending=-999;

								for (i=(int)(x[(new_n-1)*3+1]+trans_width*1.5+0.5);i>=(int)(x[(new_n-1)*3+1]+0.5); i--) /* +1 to make it the same as in matlab*/
								{
									if (i>number_of_bins)
									{
										for (i=number_of_bins;i>=(int)(x[(new_n-1)*3+1]+0.5);i--)
										{
											if ((amp_noisefree[i]>(noise)/2)&&ending<0)
											{
												ending=(int)xdata[i];
											}
										}
									}
									else
									{
										if ((amp_noisefree[i]>(noise)/2)&&ending<0)
										{
											ending=(int)xdata[i];
										}
									}
								}

								if (ending==-999)
								{
									for (i=(int)(x[(new_n-1)*3+1]+trans_width*1.5+0.5); i>=(int)(x[(new_n-1)*3+1]+0.5); i--)
									{
										if ((amp_noisefree[i]>(noise)/4)&&ending<0)
										{
											ending=(int)xdata[i];
										}
									}
								}

								for (i=1;i<=profile_bins+3*trans_width-1;i++)
								{
									if (fabs(begining_Z-height_bin[i])<(sample*vectorZ/2))
									{
										begining_bin=i;
									}
								}

								/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
						
								/*2e. AREA OF EACH BIN*/
								for (i=1; i<=begining_bin;i++)
								{
									area[i]=0;
								}
								offset_bin=begining_bin;

								
								if (offset_bin<0)
								{
									printf("\nTrouble! Pulse ID %d: profile height is %lf while waveform start is %lf \n-> increase default tree height\n",pulseID,height_bin[1],begining_Z);
									fprintf(out_sum,"\nTrouble! profile height is %lf while waveform start is %lf \n-> increase default tree height\n",height_bin[1],begining_Z);
									return(-1);
								}


								for (i=begining+1;i<=ending;i++)
								{
									offset_bin=offset_bin+1; /*adds one- pole powierzchni poprzedniego z nastepnym*/
									area[offset_bin]=((amp_noisefree[i]+amp_noisefree[i+1])/2)*(-Vector_Z);
								}

								for (i=offset_bin;i<=profile_bins+3*trans_width-2;i++)
								{
									area[i+1]=0;
								}

								for (i=1;i<=profile_bins+3*trans_width-1;i++)
								{
									energy_sum[i]=energy_sum[i]+area[i];
								}

								count_waveform_grid=count_waveform_grid+1;

								/*2f. Calulate vegetation and ground energy */
								for (i=1;i<=profile_bins+3*trans_width-1;i++)
								{
									if (height_bin[i]==0) /*find zero elevation and corresponding bin*/
									{	
										zero=i;
									}
								}

								if (last_canopy_bin_default==-999)
								{
									last_canopy=-999; 
									for (i=zero; i>=(zero-trans_width-2);i--)
									{
										if (area[i]<area[i+1]&&(area[i]<area[i-1]||(area[i]<=area[i-1]&&area[i]<area[i-2])))
										{	
											canopy_end=height_bin[i];
											last_canopy=i;
										}
									}
                
									if (last_canopy==-999)
									{
										canopy_end=height_bin[zero-trans_width/2-2];
										last_canopy=zero-trans_width/2-2;
									}
								}
								else
								{
										canopy_end=height_bin[zero-last_canopy_bin_default];
										last_canopy=zero-last_canopy_bin_default;
								}

								/*find the end of ground return in a waveform*/
								last_ground=-999; 
								for (i=zero;i<=zero+trans_width-2;i++)
								{
									if (area[i]<=area[i+1]&& area[i]<area[i-1])
									{
										ground_end=height_bin[i];
										last_ground=i;
									}
								}

								if (last_ground==-999)
								{
									ground_end=height_bin[zero+trans_width/2+2];
									last_ground=zero+trans_width/2+2;
								}

								/*calculate the energy of vegetation in a waveform*/
								energy_veg=0;
								for (i=1;i<=last_canopy;i++)
								{
									energy_veg=energy_veg+area[i];
								}

								/*calculate the energy of ground in a waveform*/
								energy_ground=0;
								for (i=last_canopy+1;i<=last_ground;i++)
								{
									energy_ground=energy_ground+area[i];
								}

								/*sum up the energy of single peak ground returns*/
								if (energy_veg==0)
								{
									mean_single_ground=mean_single_ground+energy_ground;
 									count_ground_single_grid=count_ground_single_grid+1;
								}

							} //end of waveforms falling within certain block
							else
							{
								energy_veg=-999;
								energy_ground=-999;
							} 

							NEW_N1[index]=new_n;
							energy_veg_array[index]=energy_veg;
							energy_ground_array[index]=energy_ground;
						} // end of 'if waveform is not disregarded
					}  //segment loop
				} // end of number of waveforms
            
				/*mean energy of single peak ground returns*/
				mean_single_ground_final[ygrid+1][xgrid+1]=mean_single_ground/count_ground_single_grid;
				if (mean_single_ground==0||count_ground_single_grid==0)
				{
					mean_single_ground_final[ygrid+1][xgrid+1]=-999;
				}

				
			/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
			/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

				/*For each block calculate:
				... RETURNED ENERGY
				.........reflectance ratio
				... CUMULATIVE AREA UNDER THE ENERGY GRAPH*
				... CLOSURE - NORMALIZED CUMULATIVE AREA - FROM ENERGY GRAPH
				... LEAF AREA INDEX FROM ENERGY GRAPH
				... CANOPY HEIGHT PROFILE*/

			   
				/*RETURNED ENERGY*/
				
				if (count_waveform_grid>=1)
				{
					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
							energy[j]=0;
	
					for (i=1;i<=profile_bins+3*trans_width-1;i++)
					{
						/*average bin values for the dataset*/
						energy[i]=energy_sum[i]/(double)(count_waveform_grid); /*calculate mean energy value of each column*/
					}

					if (input_YN[0]=='Y'||input_YN[0]=='y')
					{
						fprintf(out_fp1,"[%d,%d],,",ygrid+1,xgrid+1);

						for (i=1;i<=profile_bins+3*trans_width-1;i++)
							fprintf(out_fp1,"%f,",energy[i]);

						fprintf(out_fp1,"\n");
					}
        

					for (i=1;i<=profile_bins+3*trans_width-1;i++)
					{
						if (height_bin[i]==0) /*find zero elevation and corresponding bin*/
						{	
							zero=i;
						}
					}

					if (last_canopy_bin_default==-999)
					{
						last_canopy_bin=-999; /*find the end of vegetation/begining of ground return in a dataset*/
						for (i=zero; i>=(zero-trans_width-1); i--)
						{
							if (energy[i]<energy[i+1] && energy[i]<energy[i-1])
							{
								canopy_return=height_bin[i];
								last_canopy_bin=i;
							}
						}

						if (last_canopy_bin==-999)
						{
							canopy_return=height_bin[zero-trans_width/2-2];
							last_canopy_bin=zero-trans_width/2-2;
						}
					}
					else
					{
						canopy_return=height_bin[zero-last_canopy_bin_default];
						last_canopy_bin=zero-last_canopy_bin_default;
					}

					last_ground_bin=-999; /*find the end of ground return in a dataset*/
					for (i=zero; i<=profile_bins+1.5*trans_width; i++)
					{
						if (energy[i]<=energy[i+1]&& energy[i]<energy[i-1])
						{
								ground_return=height_bin[i];
								last_ground_bin=i;
						}
					}

					if (last_ground_bin==-999)
					{
						ground_return=height_bin[zero+trans_width/2+2];
						last_ground_bin=zero+trans_width/2+2;
					}

					energy_veg_total=0; /*calculate total energy of vegetation (dataset)*/
					for (i=1;i<=last_canopy_bin;i++)
					{
						energy_veg_total=energy_veg_total+energy[i];
					}

    
					energy_ground_total=0; /*calculate total energy of ground (dataset)*/
					for (i=last_canopy_bin+1;i<=last_ground_bin;i++)
					{
						energy_ground_total=energy_ground_total+energy[i];
					}

				
					/*reflectance ratio from dataset*/
					ratio_refl_total=-energy_veg_total/(energy_ground_total-mean_single_ground_final[ygrid+1][xgrid+1]);

					if (ratio_refl_total<plot_reflectance_ratio/2||(energy_ground_total-mean_single_ground_final[ygrid+1][xgrid+1]==0)||ratio_refl_total>plot_reflectance_ratio*1.5)
					{
						ratio_refl_total=plot_reflectance_ratio;
					}

					if (energy_ground_total==0||energy_veg_total==0||mean_single_ground_final[ygrid+1][xgrid+1]==-999)
					{
						ratio_refl_total=plot_reflectance_ratio;
					}
      
					/*to account for difference in reflectance*/
					for (i=last_canopy_bin+1; i<=profile_bins+3*trans_width-1;i++)
					{
						energy[i]=energy[i]*ratio_refl_total; /*to account for difference in reflectance*/
					}

					if (input_YN[0]=='Y'||input_YN[0]=='y')
						fprintf(out_fp5,"[%d,%d]: %f, %f, %f, %f, \n",ygrid+1,xgrid+1,energy_veg_total,energy_ground_total,ratio_refl_total,mean_single_ground_final[ygrid+1][xgrid+1]);
					
					refl_array[ygrid+1][xgrid+1]=ratio_refl_total;
					printf("\n Reflectance ratio in [%d,%d] is %.3f,",ygrid+1, xgrid+1,refl_array[ygrid+1][xgrid+1]);
       

					/*CUMULATIVE AREA UNDER THE ENERGY GRAPH*/

					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
							cum_area[j]=0;

					for (i=1;i<=profile_bins+3*trans_width-2;i++)
					{
						cum_area[i+1]=cum_area[i]+(vectorZ*(energy[i]+energy[i+1])/2);
					}
					total_area=cum_area[profile_bins+3*trans_width-1];

					/*CLOSURE - NORMALIZED CUMULATIVE AREA - FROM ENERGY GRAPH*/
					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
							closure[j]=0;

					total=-999;
					for (i=1;i<=profile_bins+3*trans_width-2;i++)
					{
						closure[i+1]=cum_area[i+1]/total_area;
						if (closure[i+1]==1&&total<0);
						{
							total=i+1; //a bin where the closure reaches 1
						}
					}

					/*canopy to ground energy ratio*/
					ratio=closure[last_canopy_bin]/(1-closure[last_canopy_bin]);

					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
							closure_plot[j]=0;

					for (i=1;i<=last_canopy_bin-1;i++)/*to plot*/
					{
						closure_plot[i+1]=cum_area[i+1]/total_area;
					}

					for (i=last_canopy_bin; i<=profile_bins+3*trans_width-2;i++) /*to plot*/
					{
						closure_plot[i+1]=-999;
					}

					if (input_YN[0]=='Y'||input_YN[0]=='y')
					{
						fprintf(out_fp2,"[%d,%d],,",ygrid+1,xgrid+1);
        
						for (i=1;i<=profile_bins+3*trans_width-1;i++)
							fprintf(out_fp2,"%f,",closure_plot[i]);
						fprintf(out_fp2,"\n");
					}


					/*LEAF AREA INDEX FROM ENERGY GRAPH*/
					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
							lai[j]=0;
					i=0;
					for (i=1; i<=(total-2);i++)
					{
						if (closure[i+1]!=1)
						{
							lai[i+1]=-log(1-closure[i+1]);
							total_lai=lai[i+1];
						}
						else
						{
							lai[i+1]=-999;
						}
					}

					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
							lai_plot[j]=0;

					for (i=1;i<=last_canopy_bin-1 ;i++)/*to plot*/
					{
						if (closure[i+1]!=1)
						{
							lai_plot[i+1]=-log(1-closure[i+1]);
							overall_lai=lai_plot[i+1]; /*total LAI in that grid*/
						}
						else
						{
							lai_plot[i+1]=-999;
							overall_lai=-999;
						}
					}


					LAI_array[(ygrid+1)][(xgrid+1)]=overall_lai;
					if (overall_lai!=-999)
						printf("\n LAI in [%d,%d] cell is %.3f,",ygrid+1,xgrid+1,LAI_array[ygrid+1][xgrid+1]);
					else
					{
						printf("\n LAI in [%d,%d] cell is not possible to calculate (no ground return)",ygrid+1,xgrid+1);
						fprintf(out_sum,"\n LAI in [%d,%d] cell is not possible to calculate (no ground return)",ygrid+1,xgrid+1);
						count_no_ground++; // number of cells with no ground return
					}


					for (i=last_canopy_bin;i<=profile_bins+3*trans_width-2;i++) /*to plot*/
					{
						lai_plot[i+1]=-999;
					}

					for (i=(int)total-1;i<=profile_bins+3*trans_width-2;i++)
					{
						lai[i+1]=-999;
					}

					if (input_YN[0]=='Y'||input_YN[0]=='y')
					{
						fprintf(out_fp3,"[%d,%d],,",ygrid+1,xgrid+1);
        
						for (i=1;i<=profile_bins+3*trans_width-1;i++)
							fprintf(out_fp3,"%f,",lai_plot[i]);
						fprintf(out_fp3,"\n");
					}

					/*CANOPY HEIGHT PROFILE*/
					for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
					{
							chp_grid[j]=0;
							chp_plot_grid[j]=0;
					}

					for (i=1;i<=(int)(total-2);i++);
					{
						chp_grid[i+1] = lai[i+1]-lai[i];
					}

					for (i=(int)total-1; i<=profile_bins+3*trans_width-2;i++);
					{
						chp[i+1]=-999;
						chp_grid[i+1] =-999;
					}

					for (i=1;i<=last_canopy_bin-1;i++) /*to plot*/
					{
						chp_plot_grid[i+1]=lai_plot[i+1]-lai_plot[i];
					}

					for (i=1;i<=last_canopy_bin-1;i++) /*to plot*/
					{
						if (chp_plot[i+1]==1)
						{
							chp_plot[i+1]=0;
							chp_plot_grid[i+1] = 0;
						}

					}

					for (i=last_canopy_bin;i<=profile_bins+3*trans_width-2;i++) /*to plot*/
					{
						chp_plot[i+1]=-999;
						chp_plot_grid[i+1] = -999;
					}

	

				
					for (j=1; j<=profile_bins+3*trans_width-2; j++)
					{
						if (overall_lai!=-999)
						{
							CHP_array[((ygrid)*xgrid_total+(xgrid+1))][j]=chp_plot_grid[j];
						}
						else
						{
							CHP_array[((ygrid)*xgrid_total+(xgrid+1))][j]=0;
						}
					}

					for (j=1; j<=profile_bins+3*trans_width-2; j++)
					{
						if (CHP_array[((ygrid)*xgrid_total+(xgrid+1))][j]==-999)
						{
							CHP_array[((ygrid)*xgrid_total+(xgrid+1))][j]=0;
						}
					}

					if (input_YN[0]=='Y'||input_YN[0]=='y')
					{
						fprintf(out_fp4,"[%d,%d],,",ygrid+1,xgrid+1);	
						for (j=1; j<=profile_bins+3*trans_width-2; j++)
						{
							fprintf(out_fp4,"%f,",CHP_array[((ygrid)*xgrid_total+(xgrid+1))][j]);
						}
						fprintf(out_fp4,"\n");
					}
					printf("\n");
				}
				else
				{
					LAI_array[(ygrid+1)][(xgrid+1)]=0;
					for (j=1; j<=profile_bins+3*trans_width-2; j++)
					{
							CHP_array[((ygrid)*xgrid_total+(xgrid+1))][j]=0;
					}

					count_no_data++;
					printf("\n WARNING: No data in cell [%d,%d]\n",ygrid+1,xgrid+1);

					if (input_YN[0]=='Y'||input_YN[0]=='y')
					{
						fprintf(out_fp1,"[%d,%d],, no data\n",ygrid+1,xgrid+1);
						fprintf(out_fp2,"[%d,%d],, no data\n",ygrid+1,xgrid+1);
						fprintf(out_fp3,"[%d,%d],, no data\n",ygrid+1,xgrid+1);
						fprintf(out_fp4,"[%d,%d],, no data\n",ygrid+1,xgrid+1);
						fprintf(out_fp5,"[%d,%d]: No_data, No_data, No_data, No_data, \n",ygrid+1,xgrid+1);
					}

					fprintf(out_sum,"[%d,%d],, WARNING: no data in the cell\n",ygrid+1,xgrid+1);
				}
			}
		}//end of gridcell


		//////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////
		
		if (count_no_data==ygrid_total*ygrid_total)
		{
			printf("\nNo data in the specified area");
			fprintf(out_sum,"\nNo data in the specified area");
 			return(-1);
		}

		// find a maximum LAI grid value in the dataset
		LAI_maximum=0;
		for (i=1;i<=ygrid_total;i++)
		{
			for (j=1;j<=xgrid_total;j++)
			{
				if (LAI_maximum<LAI_array[i][j])
				{
					LAI_maximum=LAI_array[i][j];
				}
			}
		}

		for (i=1;i<=ygrid_total;i++)
		{
			for (j=1;j<=xgrid_total;j++)
			{
				if (LAI_array[i][j]!=-999)
				{}
				else
				{
					LAI_array[i][j]=LAI_maximum;
				}
			}
		} 

		fprintf(out_fp8,"%lf, %lf\n", xcell_size, ycell_size);
		fprintf(out_fp8,"%lf, %lf\n", xdim, ydim);
		fprintf(out_fp8,"%lf, %lf\n", x_ll, y_ll);
		printf("\n");
		
		/*Output LAI map values*/
		for (i=1;i<=ygrid_total;i++)
		{
			printf("LAI row %d,", i);
			for (j=1;j<=xgrid_total;j++)
			{
				fprintf(out_fp8,"%f,", LAI_array[i][j]);
				
				printf("%f,", LAI_array[i][j]);
			}
			fprintf(out_fp8,"\n");
			printf("\n");
		}

		/*Output Reflectance values*/
		if (input_YN[0]=='Y'||input_YN[0]=='y')
			fprintf(out_fp5,"\n Reflectance array\n");
		printf("\n");
		for (i=1;i<=ygrid_total;i++)
		{
			printf("Row_%d,", i);

			if (input_YN[0]=='Y'||input_YN[0]=='y')
				fprintf(out_fp5,"Row_%d,", i);

			for (j=1;j<=xgrid_total;j++)
			{
				if (input_YN[0]=='Y'||input_YN[0]=='y')
					fprintf(out_fp5,"%f,", refl_array[i][j]);

				printf("%f,", refl_array[i][j]);
			}
			if (input_YN[0]=='Y'||input_YN[0]=='y')
				fprintf(out_fp5,"\n");
			printf("\n");
		}

		/*calculate total LAI*/
		Total_LAI=0;
		for (i=1;i<=ygrid_total;i++)
		{
			for (j=1;j<=xgrid_total;j++)
			{
				Total_LAI=Total_LAI+LAI_array[i][j];	
			}
		}
		if ((ygrid_total*xgrid_total)-count_no_data!=0)
		{
			Total_LAI=Total_LAI/((ygrid_total*xgrid_total)-count_no_data);
			fprintf(out_sum,"\n\nTotal site LAI (LAI1): \n\tsaturated cells replaced with \n\tthe maximum dataset LAI value: \t%f",Total_LAI);
			
			if (input_YN[0]=='Y'||input_YN[0]=='y')
				fprintf(out_fp3,"\n\nTotal site LAI (LAI1): \n\tsaturated cells replaced with \n\tthe maximum dataset LAI value: \t%f",Total_LAI);
			
			printf("\nTotal site LAI %f,", Total_LAI);
		}
		else
		{
			fprintf(out_sum,"\n\n WARNING: No data in the selected area!!!");
			fclose(out_sum);
			return(-1);
		}
		
		//alternative LAI to avoid too hight LAI values 
		//(assigning to the cells with no ground return a mean value of total LAI calculated for the extreem case (TOTAL_LAI) - with the maximum value)
		for (i=1;i<=ygrid_total;i++)
		{
			for (j=1;j<=xgrid_total;j++)
			{
				if (LAI_maximum==LAI_array[i][j] &&count_no_ground>0)
				{
					LAI_array_corrected[i][j]=Total_LAI;
				}
				else
				{
					LAI_array_corrected[i][j]=LAI_array[i][j];
				}
			}
		}

		Total_LAI_corrected=0;
		for (i=1;i<=ygrid_total;i++)
		{
			for (j=1;j<=xgrid_total;j++)
			{
				Total_LAI_corrected=Total_LAI_corrected+LAI_array_corrected[i][j];	
			}
		}

		Total_LAI_corrected=Total_LAI_corrected/((ygrid_total*xgrid_total)-count_no_data);
		fprintf(out_sum,"\n\nTotal site LAI (LAI2): \n\tsaturated cells replaced with \n\tthe mean total site LAI1 value: %f",Total_LAI_corrected);
		
		if (input_YN[0]=='Y'||input_YN[0]=='y')
			fprintf(out_fp3,"\n\nTotal site LAI (LAI2): \n\tsaturated cells replaced with \n\tthe mean total site LAI1 value: %f",Total_LAI_corrected);
		
		printf("\nTotal site LAI %f,", Total_LAI_corrected);

		
		/*calculate total reflectance ratio*/
		Total_refl=0;
		for (i=1;i<=ygrid_total;i++)
		{
			for (j=1;j<=xgrid_total;j++)
			{
				Total_refl=Total_refl+refl_array[i][j];	

			}
		}
		Total_refl=Total_refl/((ygrid_total*xgrid_total)-count_no_data);
		
		if (input_YN[0]=='Y'||input_YN[0]=='y')
			fprintf(out_fp5,"\n\nTotal site reflectance ratio %f",Total_refl);

		printf("\nTotal site reflectance ratio %f,", Total_refl);

		////////////////////////////////////////////////////////////////
		//CHP
		for (j=0;j<=profile_bins+3*trans_width-1;j++) /*zero*/
		{
			chp_total_sum[j]=0;
			chp_total[j]=0;		
			chp_total_norm[j]=0;
		}
		chp_total_sum_sum=0;

		//add up the CHP values bin-wise (normalization after adding up all cells)
		for (j=1; j<=profile_bins+3*trans_width-1; j++)
		{
			for (i=1; i<=ygrid_total*xgrid_total; i++)
			{
				chp_total_sum[j]=chp_total_sum[j] + CHP_array[i][j];
			}
		}

		
		//calculate average bin value - unnormalized CHP
		for (j=1; j<=profile_bins+3*trans_width-1; j++)
			chp_total[j]=chp_total_sum[j]/(ygrid_total*xgrid_total-count_no_ground-count_no_data);

		
		// calculate total CHP value (for normalization) = total LAI
		for (j=1; j<=profile_bins+3*trans_width-1; j++)
			chp_total_sum_sum=chp_total_sum_sum+chp_total[j];
		
		//calculate final normalized CHP
		for (j=1; j<=profile_bins+3*trans_width-1; j++)
			chp_total_norm[j]=chp_total[j]/chp_total_sum_sum;


		for(j=zero; j<=profile_bins+3*trans_width-1; j++)
		{
			chp_total[j]=-999;
		}

		
		fprintf(out_sum,"\n\nTotal site LAI (LAI3): \n\twith saturated cells removed: \t%f", chp_total_sum_sum);
		fprintf(out_sum,"\n\nDifference between LAI1 and LAI3: %f\nDifference between LAI2 and LAI3: %f", Total_LAI-chp_total_sum_sum, Total_LAI_corrected-chp_total_sum_sum);

		if (input_YN[0]=='Y'||input_YN[0]=='y')
		{
			fprintf(out_fp3,"\n\nTotal site LAI (LAI3): \n\twith saturated cells removed: \t%f", chp_total_sum_sum);
			fprintf(out_fp3,"\n\nDifference between LAI1 and LAI3: %f\nDifference between LAI2 and LAI3: %f", Total_LAI-chp_total_sum_sum, Total_LAI_corrected-chp_total_sum_sum);
		}
		
		fprintf(out_sum,"\n\nNumber of cells with no ground return: %d", count_no_ground);
		fprintf(out_sum,"\n\nNumber of cells with no data: %d", count_no_data);

		fprintf(out_sum,"\n\nTotal site reflectance ratio %f",Total_refl);
		//fprintf(out_fp4,"\n\nTotal CHP,,,");

		for(j=zero; j<=profile_bins+3*trans_width-1; j++)
		{
			chp_total_norm[j]=-999;
		}

		fprintf(out_fp9,"%d\n",profile_bins+3*trans_width-1);

		for(j=1; j<=profile_bins+3*trans_width-1; j++)
		{
			fprintf(out_fp9,"%f,", chp_total_norm[j]);
		}

	}

	fprintf(out_sum,"\n\n*************************************************************************");
	fprintf(out_sum,"\n\nPROCESSING COMPLETED!");

	printf("\n\nPROCESSING COMPLETED!");

	if (input_YN[0]=='Y'||input_YN[0]=='y')
	{
		fclose(out_fp1);
		fclose(out_fp2);
		fclose(out_fp3);
		fclose(out_fp4);
		fclose(out_fp5);
	}

	fclose(out_sum);
	fclose(out_fp8);
	fclose(out_fp9);
	fclose(in_fp);
	fclose(in_fp2);
	fclose(in_fp3);
	fclose(in_fp4);
	
	return(0);	
}
	
