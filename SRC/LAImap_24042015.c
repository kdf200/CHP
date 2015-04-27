//***********************************************************************
//**               Module 4a of Canopy Height Program                  **
//***********************************************************************
//** Produces LAIe map using PLplot library                            **
//**                                                                   **
//**                                                                   **
//** Authors: Karolina Fieber, Ian Davenport                           **
//** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **
//** University of Reading, UK              November 2014 - March 2015 **
//***********************************************************************

#include "plcdemos.h"

// Fundamental settings.  See notes[] for more info.
static int ns      = 20;        // Default number of shade levels

// Function prototypes

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

static void
f2mnmx( double **f, int xgrid_total, int ygrid_total, double *fmin, double *fmax );
// For now, don't show the colorbars while we are working out the API.


static int colorbar = 1;


static double  x_LAI[] = { 0.1, 0.5, 0.5, 0.1 };
static double  y_LAI[] = { 0.1, 0.1, 0.5, 0.5 };

//--------------------------------------------------------------------------
// main
//
// Does several shade plots using different coordinate mappings.
//--------------------------------------------------------------------------

int
main( int argc, const char *argv[] )
{
    int        i, j;
    double      *clevel, *shedge;

    double      fill_width = 2., cont_width = 0.;
    double      colorbar_width, colorbar_height;
    int      cont_color = 0;
#define NUM_AXES    1
    int      n_axis_opts  = NUM_AXES;
    const char *axis_opts[] = {
        "bcvtm",
    };
    int      num_values[NUM_AXES];
    double      *values[NUM_AXES];
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
        "LAIe",
    };

	FILE	*in_fp;
	double 	**LAI_array_corrected;
	int	ygrid_total,xgrid_total;
	double	xcell_size,ycell_size,x_ll,y_ll,xdim,ydim;

	double	LAI_maximum, LAI_minimum;
	double	empty;
	
	double	reminder,reminderY;
	char	identifier_outX[25], identifier_outY[25],characX[50], characY[50];
	char	in_filename[200],filename[200];
	char         string[30];

	printf("***********************************************************************\n");
	printf("**               Module 4a of Canopy Height Program                  **\n");
	printf("***********************************************************************\n");
	printf("** Produces LAIe map using PLplot library                            **\n");
	printf("**                                                                   **\n");
	printf("**                                                                   **\n");
	printf("** Authors: Karolina Fieber, Ian Davenport                           **\n");
	printf("** k.fieber@reading.ac.uk, i.j.davenport@reading.ac.uk               **\n");
	printf("** University of Reading, UK                              March 2015 **\n");
	printf("***********************************************************************\n\n");


// Load colour palettes
    plspal0( "cmap0_black_on_white.pal" );
    plspal1( "cmap1_gray.pal", 1 );

// Initialize plplot
	plsdev("wingcc");
    

// Allocate data structures

    clevel = (double *) calloc( (size_t) ns, sizeof ( double ) );
    shedge = (double *) calloc( (size_t) ( ns + 1 ), sizeof ( double ) );

// Set up data array
	printf("Enter LAI input filename (full without extension): ");
	scanf("%s",&filename);
	
	sprintf(in_filename,"%s.txt",filename);
	printf("\n\tReading from file <<%s>>\n",in_filename);
	in_fp=fopen(in_filename,"r");
	if (!in_fp)
	{
		printf("\tCan't open input data file <%s>\n",in_filename);
 		return(-1);
	}

	xcell_size=readindouble(in_fp);
	ycell_size=readindouble(in_fp);
	xdim=readindouble(in_fp);
	ydim=readindouble(in_fp);
	x_ll=readindouble(in_fp);
	y_ll=readindouble(in_fp);
	
	ygrid_total=(int)(ydim/ycell_size); /*number of grid cell -rows*/
	xgrid_total=(int)(xdim/xcell_size); /*number of grid cell -columns*/

	LAI_array_corrected=(double **)calloc(ygrid_total+2,sizeof(double));
	for (j=0;j<=ygrid_total+1;j++)
		LAI_array_corrected[j]=(double *)calloc(xgrid_total+1,sizeof(double));

	for (i=0;i<=ygrid_total+1;i++)
		for (j=0;j<=xgrid_total+1;j++)
		  {
					LAI_array_corrected[i][j]=(double)(0);
		  }

	LAI_maximum=-999;

	//Read in LAI array
	for (i=0;i<=ygrid_total-1;i++)
	{			
		for (j=0;j<=xgrid_total-1;j++)
		{
			LAI_array_corrected[i][j]=readindouble(in_fp);
		}
		empty=readindouble(in_fp);
	}


    f2mnmx( LAI_array_corrected, xgrid_total-1, ygrid_total-1, &LAI_minimum, &LAI_maximum );
	LAI_minimum=0;
    for ( i = 0; i < ns; i++ )
        clevel[i] = LAI_minimum + ( LAI_maximum - LAI_minimum ) * ( i + 0.5 ) / (double) ns;

    for ( i = 0; i < ns + 1; i++ )
        shedge[i] = LAI_minimum + ( LAI_maximum  - LAI_minimum ) * (double) i / (double) ns;


// Plot using 2d coordinate transform
	plinit();
// Load colour palettes

    plspal0( "cmap0_black_on_white.pal" );
    plspal1( "cmap1_blue_red.pal", 1 );

	plschr( 0.,0.7);
	plenv(x_ll,x_ll+xdim,y_ll,y_ll+ydim,1,0);

	// Set up coordinate grids
	for (i=0;i<=ygrid_total-1;i++)
	{
		y_LAI[0]=y_ll+(i)*ycell_size;
		y_LAI[1]=y_ll+(i)*ycell_size;
		y_LAI[2]=y_ll+(i+1)*ycell_size;
		y_LAI[3]=y_ll+(i+1)*ycell_size;

		for (j=0;j<=xgrid_total-1;j++)
		{					
			x_LAI[0]=x_ll+(j)*xcell_size;
			x_LAI[1]=x_ll+(j+1)*xcell_size;
			x_LAI[2]=x_ll+(j+1)*xcell_size;
			x_LAI[3]=x_ll+(j)*xcell_size;

			//Assign colour to every cell according to LAI value
			plcol1(LAI_array_corrected[i][j]/ (LAI_maximum));	
			plfill( 4, x_LAI, y_LAI);
		}
		//printf("\n");
	}

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

    plcol0( 2 );

	plschr( 0.,1.0 );
	plmtex("b",2.5, 0.5, 0.5, "Easting [m]");
	plmtex("l",4.8, 0.5, 0.5, "Northing [m]");
	plmtex("t",2.0, 0.5, 0.5, "LAI MAP ");

	plfont( 1 );
    plschr( 0., 0.7 );
	sprintf(string, "file: %s.txt", filename);
	plmtex("r", 8, 0.5, 0.5, string);
	sprintf(string, "grid: %.1lfm by %.1lfm", xcell_size, ycell_size);
	plmtex("t", 4.0, 1.0, 1.0, string);
	sprintf(string, "area: %.1lfm by %.1lfm", xdim, ydim );
	plmtex("t", 2.5, 1.0, 1.0, string);
	sprintf(string, "LLC: %.0lfm, %.0lfm", x_ll, y_ll );
	plmtex("t", 1.0, 1.0, 1.0, string);


// Clean up

    free( (void *) clevel );
    free( (void *) shedge );

    plend();

    exit( 0 );
}

//--------------------------------------------------------------------------
// f2mnmx
//
// Returns min & max of input 2d LAIe array.
//--------------------------------------------------------------------------

static void
f2mnmx( double **f, int nnx, int nny, double *fnmin, double *fnmax )
{
    int i, j;

    *fnmax = f[0][0];
    *fnmin = *fnmax;

    for ( i = 0; i < nny+1; i++ )
    {
        for ( j = 0; j < nnx+1; j++ )
        {
            *fnmax = MAX( *fnmax, f[i][j] );
            *fnmin = MIN( *fnmin, f[i][j] );
        }
    }
}
