/*
gcc mandelbrot_omp.c -o mandelbrot_omp -lglut -lGLU -lGL -fopenmp
*/
#include <stdlib.h>
#include <GL/glut.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

/* Defaut data via command line */
/* Can enter other values via command line arguments */

#define CENTERX -0.5
#define CENTERY 0.5
#define HEIGHT 0.5
#define WIDTH 0.5
#define MAX_ITER 100

/* N x M array to be generated */

#define N 500
#define M 500

/* Number of threads*/

#define NUM_THREADS 8

float height = HEIGHT; /* size of window in complex plane */
float width = WIDTH;
float cx = CENTERX; /* center of window in complex plane */
float cy = CENTERY;
int max = MAX_ITER; /* number of interations per point */

int n=N;
int m=M;

/* Use unsigned bytes for image */

GLubyte image[N][M];

/* Complex data type and complex add, mult, and magnitude functions */
/* Probably not worth overhead */

typedef float complex[2];

void add(complex a, complex b, complex p)
{
    p[0]=a[0]+b[0];
    p[1]=a[1]+b[1];
}

void mult(complex a, complex b, complex p)
{
    p[0]=a[0]*b[0]-a[1]*b[1];
    p[1]=a[0]*b[1]+a[1]*b[0];
}

float mag2(complex a)
{
    return(a[0]*a[0]+a[1]*a[1]);
}

void form(float a, float b, complex p)
{
    p[0]=a;
    p[1]=b;
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(n,m,GL_COLOR_INDEX, GL_UNSIGNED_BYTE, image);
    glFlush();
}


void myReshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (w <= h)
        gluOrtho2D(0.0, 0.0, (GLfloat) n, (GLfloat) m* (GLfloat) h / (GLfloat) w);
    else
        gluOrtho2D(0.0, 0.0, (GLfloat) n * (GLfloat) w / (GLfloat) h,(GLfloat) m);
    glMatrixMode(GL_MODELVIEW);
    display();
}

void myinit()
{
    float redmap[256], greenmap[256],bluemap[256];
    int i;

    glClearColor (1.0, 1.0, 1.0, 1.0);
    gluOrtho2D(0.0, 0.0, (GLfloat) n, (GLfloat) m);

    /* Define pseudocolor maps, ramps for red and blue,
       random for green */


    for(i=0; i<256; i++)
    {
        redmap[i]=i/255.;
        greenmap[i]=drand48();
        bluemap[i]=1.0-i/255.;
    }

    glPixelMapfv(GL_PIXEL_MAP_I_TO_R, 256, redmap);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_G, 256, greenmap);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_B, 256, bluemap);
}


main(int argc, char *argv[])
{
    int i, j, k;
    float x, y, v;
    complex c0, c, d;
    double elapsedTime;
    struct timeval start, end;

    if(argc>1) cx = atof(argv[1]); /* center x */
    if(argc>2) cy = atof(argv[2]);  /* center y */
    if(argc>3) height=width=atof(argv[3]); /* rectangle height and width */
    if(argc>4) max=atoi(argv[4]); /* maximum iterations */

    omp_set_num_threads(NUM_THREADS);
    gettimeofday(&start, NULL);
    #pragma omp parallel for default(shared) private(x, y, c, c0, v, j, d, k)   //Because these items are declared outside so we must use private explicitly
    for (i=0; i<n; i++)
        for(j=0; j<m; j++)
        {

            /* starting point */

            x= i *(width/(n-1)) + cx -width/2;
            y= j *(height/(m-1)) + cy -height/2;

            form(0,0,c);
            form(x,y,c0);

            /* complex iteration */

            for(k=0; k<max; k++)
            {
                mult(c,c,d);
                add(d,c0,c);
                v=mag2(c);
                if(v>4.0) break; /* assume not in set if mag > 4 */
            }

            /* assign gray level to point based on its magnitude */
            if(v>1.0) v=1.0; /* clamp if > 1 */
            image[i][j]=255*v;
        }
    gettimeofday(&end, NULL);
    elapsedTime = (end.tv_sec - start.tv_sec) * 1000;
    elapsedTime = (end.tv_usec - start.tv_usec) / 1000;
    printf("OMP_Time: %f ms\n", elapsedTime);


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(N, M);
    glutCreateWindow("mandlebrot");
    myinit();
    glutReshapeFunc(myReshape);
    glutDisplayFunc(display);

    glutMainLoop();


}
