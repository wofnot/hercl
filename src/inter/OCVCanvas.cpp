#include "OCVCanvas.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <iostream>
#include <memory>

// currently is square rather than toroidal canvas.
// toroidal was nice but drawing lines across borders is ~tricky

using namespace std;
using namespace cv;


#define SCALE 1.0

#define LINE_TYPE 8

void drawLine(Canvas canvas, double xa, double ya, double xb, double yb);

double clamp(double l, double x, double u);

float s;

Canvas newCanvas(int width, int height, bool colour) {
  Canvas canvas = make_shared<canvasStruct>();
  
  canvas->width = width;
  canvas->height = height;
  
  
  canvas->colour = colour;
  
  //; // 8-bit greyscale

  if (canvas->colour) {
   canvas->pixels = Mat(SCALE*canvas->height, SCALE*canvas->width, CV_8UC3);
  } else { 
   canvas->pixels = Mat(SCALE*canvas->height, SCALE*canvas->width, CV_8UC1);
  }
  
  resetCanvas(canvas);
  
  return canvas;
}

Canvas loadPGM(FILE *fi) {
  char line[100];
  int line_size = 100;
  int w,h;
  
  Canvas canvas = (Canvas)0;
  
  int i = 0; // pixel index
  int c,p;
  
  int n = 0;
  
  while (fgets(line, line_size, fi) != NULL)  {
    if (line[0] == '#') continue;
    if (n == 0) {
      assert(line[0] == 'P');
      assert(line[1] == '2');
    } else if (n == 1) {
      // read w,h
      c = sscanf(line, "%d %d", &w, &h);
      assert(c == 2);
      canvas = newCanvas(w,h);
    } else if (n == 2) {
      // read max value (ignore)
      c = sscanf(line, "%d", &p);
      assert(c == 1);
      assert(p == 255);
    } else {
      // read pixel
      c = sscanf(line, "%d", &p);
      assert(c == 1);
      assert(p < 256);
      assert(i < w*h);
      canvas->pixels.at<unsigned char>(i/w,i%w) = p;
      i++;
    }
    n++;
  }
  assert(canvas != (Canvas)NULL);
  return canvas;
}

double t(double x) {
    return (x+0.5)*SCALE;
}

void savePNG(Canvas canvas, char *fn) {
  //cout << "savePNG...\n";
  if (canvas->colour) {
    //cout << "in colour\n";
    Mat bgr;
    cvtColor(canvas->pixels, bgr, CV_HLS2BGR);
    //imshow("hi", canvas->pixels);
   // waitKey(0);
    imwrite(string(fn), bgr);
  } else {
    imwrite(string(fn), canvas->pixels);
  }
}

Mat getMatrix(Canvas canvas) {
  return canvas->pixels;
}
int getWidth(Canvas canvas) {
  return canvas->width;
} 
int getHeight(Canvas canvas) {
  return canvas->height;
}
int isPenDown(Canvas canvas) {
  return canvas->isPenDown;
}
void getPenPos(Canvas canvas, double *x, double *y) {
  *x = canvas->penX;
  *y = canvas->penY;
}

void getPenSize(Canvas canvas, int *size) {
  *size = canvas->penSize;
}
void getPenColour(Canvas canvas, colour_t *colour) {
  *colour = canvas->penColour;
}

void setPenSize(Canvas canvas, int size) {
  assert(size >= 0);
  if (size == 0) size = 1;
  canvas->penSize = size;
}

void togglePen(Canvas canvas) {
  canvas->isPenDown = 1 - canvas->isPenDown;
  assert(canvas->isPenDown == 1 || canvas->isPenDown == 0);
  if (canvas->isPenDown) {
    movePen(canvas, 0, 0);
  }
}

void setPenColour(Canvas canvas, colour_t colour) {
  canvas->penColour = colour;
}

void setPenPos(Canvas canvas, double x, double y) {
  double oldX = canvas->penX;
  double oldY = canvas->penY;
  canvas->penX = clamp(0, x ,canvas->width - 1);
  canvas->penY = clamp(0, y ,canvas->height - 1);
  
  if (canvas->isPenDown) {
    drawLine(canvas, oldX, oldY, canvas->penX, canvas->penY);
  }
}

void movePen(Canvas canvas, double dx, double dy) {
  setPenPos(canvas, canvas->penX+dx, canvas->penY+dy);
}

void resetCanvas(Canvas canvas) {
  if (canvas->colour) {
    canvas->pixels = Scalar(0,255,255); // white
    canvas->penColour = Scalar(0,128,255); // red ?
  } else {
    canvas->penColour = Scalar(0);
    canvas->pixels = Scalar(255);
  }
  
  canvas->isPenDown = 1;
  canvas->penX = canvas->width/2;
  canvas->penY = canvas->height/2;
  canvas->penSize = 3;

}

void showCanvas(Canvas canvas) {
  int x,y;
  for (y = 0; y < canvas->height; y++) {
    for (x = 0; x < canvas->width; x++) {
      printf(canvas->pixels.at<unsigned char>(y,x) ? "##" : "..");
    }
    printf("\n");
  }
}

void freeCanvas(Canvas canvas) {
  //free(canvas);
}


void drawLine(Canvas canvas, double x0, double y0, double x1, double y1) {
  line(canvas->pixels, Point(t(x0),t(y0)), Point(t(x1),t(y1)), canvas->penColour, SCALE*canvas->penSize, LINE_TYPE);
}

void drawCircle(Canvas canvas, double x, double y, double r) {
  circle(canvas->pixels, Point(t(x),t(y)), SCALE*r, canvas->penColour, -1, LINE_TYPE);
}

int mod (int x, int N) {
  int q = x%N;
  if (q < 0) q+=N;
  return q;
}

double clamp(double l, double x, double u) {
  if (x < l) return l;
  else if (x > u) return u;
  else return x;
}
