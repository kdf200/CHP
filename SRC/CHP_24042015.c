//***********************************************************************
//**               Module 4b of Canopy Height Program                  **
//***********************************************************************
//** Produces CHP plot using PLplot library                           **
//**                                                                   **
//**                                                                   **
//** Authors: Karolina Fieber, Ian Davenport                           **
//** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **
//** University of Reading, UK              November 2014 - March 2015 **
//***********************************************************************

#include "plcdemos.h"

// Function prototypes

double readindouble(FILE *fp)
{
	double scanned_number;
	char instring[200];
	int i;
	instring[0]=' ';
	for (i=1;(instring[i-1]!=',')&&(instring[i-1]!='\n')&&(i<99);i++)
		{
			instring[i]=fgetc(fp);
			//printf("%c",instring[i]);
	}
	if (i<99)
	{
		instring[i-1]='\0';
		sscanf(instring,"%lf",&scanned_number);
	}
	else 
		scanned_number=-999;
	return scanned_number;
}

//--------------------------------------------------------------------------
// main
//
// Does CHP plot
//--------------------------------------------------------------------------

int
main( int argc, const char *argv[] )
{
    int     i;

	FILE	*in_fp;

	double	height_bin[5001], CHP_total[5001];
	double	chp_max,xmax,ymax, xmin,vmin,vmax,profile_length, empty;
	
	char	in_filename[200],filename[200],string[30];
	int		count;

	printf("***********************************************************************\n");
	printf("**               Module 4b of Canopy Height Program                  **\n");
	printf("***********************************************************************\n");
	printf("** Produces CHP plot using PLplot library                            **\n");
	printf("**                                                                   **\n");
	printf("**                                                                   **\n");
	printf("** Authors: Karolina Fieber, Ian Davenport                           **\n");
	printf("** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	printf("** University of Reading, UK                              March 2015 **\n");
	printf("***********************************************************************\n\n");

// Load colour palettes
    plspal0( "cmap0_black_on_white.pal" );

// specify the device
	plsdev("wingcc");
// Set up data array
	printf("\nEnter CHP file name (exluding extension): ");
	scanf("%s",&filename);


	// open the file for reading
	sprintf(in_filename,"%s.txt", filename);
	in_fp=fopen(in_filename,"r");
	if (!in_fp)
	{
		printf("\n\tCan't open input data file <%s>\n",in_filename);
 		return(-1);
	}
	printf("\n\tReading from file <<%s>>\n",in_filename);

	// initialize the arrays
	for (i=0; i<5000; i++)
	{
		CHP_total[i]=0;
		height_bin[i]=0;
	}

	profile_length=readindouble(in_fp);

	//Read in Height bin profile
	for (i=0;i<= (int)profile_length-1;i++)
	{	
		height_bin[i]=readindouble(in_fp);
		//if (height_bin[i]==(-952))
			//break;
	}
	
	count=0;
	chp_max=0;
	empty=readindouble(in_fp);

	profile_length=readindouble(in_fp);

	//Read in CHP profile
	for (i=0;i<=(int)profile_length-1;i++)
	{	
		CHP_total[i]=readindouble(in_fp);
		count++;

		if(chp_max<CHP_total[i])
			chp_max=CHP_total[i];

		if (CHP_total[i]==-999)
		{
			CHP_total[i]=0;	
			//break;
		}
	}

	xmin=0.0;
	xmax=(chp_max+chp_max*0.1);

	ymax=(double) height_bin[0]+1;
	
	// Plot using 2d coordinate transform
	plinit();
	plcol0( 1 );
	plenv( (double)0.0, (double) xmax,(double) 0.0, (double) ymax, 0, 0 );


	plcol0( 3 );
	pllab( "Fraction of Foliage [%]", "Elevation [m]", "Canopy Height Profile" );

	plfont( 1 );
    plschr( 0., 0.7 );
	sprintf(string, "%s.txt", filename);
	plmtex("r", 5, 0.5, 0.5, string);

	plline( (count+2),CHP_total, height_bin);
	////////////////////////////////////////////////////
	printf("\n\tPlotting from file <<%s>>\n",in_filename);
    plend();
    exit( 0 );
}





