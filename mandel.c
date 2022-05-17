
#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void * compute_image(void * arg);

//parameters for compute image which will be used in thread creation
struct parameters
{
    struct bitmap *bm;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    int max;
    int id;
    int threads;
};

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
    printf("-n          The amount of threads used to run the program. (default=1)\n");
    printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}
int main( int argc, char *argv[] )
{
    struct timeval begin_time;
    struct timeval end_time;
    gettimeofday(&begin_time,NULL);
    char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;

    //amound of threads to be built. 1 thread if no input
    int threads = 1;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:n:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
            case 'n'://Case statement used to get number of threads inputed by user
                threads = atoi(optarg);
                break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

    //pthread creation. Will split up threads in terms of the rows of the image.
    pthread_t thrarr[threads];
    struct parameters param[threads];
	for(int x =0; x< threads; x++)
    {
        param[x].bm = bm;
        param[x].xmin = xcenter-scale;
        param[x].xmax = xcenter+scale;
        param[x].ymin = ycenter-scale;
        param[x].ymax = ycenter+scale;
        param[x].max = max;
        param[x].id = x;
        param[x].threads = threads;
        pthread_create(&thrarr[x],NULL,compute_image,(void *)&param[x]);
    }
    for(int x =0;x < threads;x++)
    {
        pthread_join(thrarr[x],NULL);
    }
    // Compute the Mandelbrot image
	//compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}
    gettimeofday(&end_time,NULL);
    int time_to_execute = ( end_time.tv_sec * 1000000 + end_time.tv_usec ) -
                         ( begin_time.tv_sec * 1000000 + begin_time.tv_usec );

  printf("This code took %d microseconds to execute\n", time_to_execute);
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/
//Added id and threads in order to know which rows to make in the image
void * compute_image(void * arg)
{
	int i,j;
    //parameters from pointer
    struct parameters * params = (struct parameters *) arg;
    //setting variables from parameters
    struct bitmap *bm = params->bm;
    double xmin = params->xmin;
    double xmax = params->xmax;
    double ymin = params->ymin;
    double ymax = params->ymax;
    int max = params->max;
    int threads = params->threads;
    int id = params->id;

    int width = bitmap_width(bm);
    int height = bitmap_height(bm);
    //used to see rows needed for individual threads
    int heightadj = height/threads;


	// For every pixel in the image...
	for(j=id*heightadj;j<heightadj*(id+1);j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
    return NULL;
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}




