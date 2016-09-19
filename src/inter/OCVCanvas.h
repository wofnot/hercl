#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <memory>

#define COLOUR

using namespace cv;
using namespace std;


typedef Scalar colour_t ; // if colour: (h,l,s), otherwise (l)


struct canvasStruct {
  Mat pixels;
  int width;
  int height;
  
  bool colour;

  int isPenDown;
  double penX;
  double penY;
  int penSize;
  colour_t penColour;
};

typedef shared_ptr<canvasStruct> Canvas;

Canvas newCanvas(int width, int height, bool colour = false);

Canvas loadPGM(FILE *fi);

void savePNG(Canvas canvas, char *fn);

Mat getMatrix(Canvas canvas);
int getWidth(Canvas canvas);
int getHeight(Canvas canvas);

int isPenDown(Canvas canvas);
void getPenPos(Canvas canvas, double *x, double *y);
void getPenSize(Canvas canvas, int *size);
void getPenColour(Canvas canvas, colour_t *colour);

void setPenSize(Canvas canvas, int size);

void togglePen(Canvas canvas);

void movePen(Canvas canvas, double dx, double dy);


void setPenColour(Canvas canvas, colour_t colour);

void setPenPos(Canvas canvas, double x, double y);


void drawCircle(Canvas canvas, int x, int y, int r);

void resetCanvas(Canvas canvas);

void showCanvas(Canvas canvas);

void freeCanvas(Canvas canvas);
