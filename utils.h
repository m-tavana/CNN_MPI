#ifndef UTILS_H
#define UTILS_H

//the user desired level of parallelism
#define UserParallelism 32

void read_weight1(const char filename[], int size, float matrix[]);
void read_bias1(const char filename[], int length, float vector[]);
void read_image_pgm(unsigned char image[], char filename[],
                    int imageWidth, int imageHeight);
void write_image_pgm(unsigned char image[], const char filename[], int imageWidth, int imageHeight);

#endif//UTILS_H
