#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_BUFFER_SIZE 100  //Line buffer size for read write 

/*******************************************************************************
 * Input   : char array containing the filename containing all weights, number of weights
 * Output  : array matrix filled with weights for each feature map
 * Procedure: read all weights from file and strore in array
 ******************************************************************************/
void read_weight1(const char filename[], int size, float matrix[]) {
  FILE* finput;
    
  finput = fopen(filename , "rb" );
  if (finput==NULL) {fputs ("File error",stderr); exit (13);}
  
  fread(matrix, sizeof(float), size, finput);
  fclose(finput);
}

/************************************************************************************
 * Function: void read_bias(char filename[], int length, float vector[])
 * Input   : char array containing the filename and location for reading, number of bias values this
                is the same as the number of output featuremaps, pointer for output
                * Output  : vector filled with bias weights for each feature map
                * Procedure: read bias weights from file normalize to uchar range and strore on correct possition
                ************************************************************************************/
void read_bias1(const char filename[], int length, float vector[]) {
  int i;
  FILE* finput;
  
  finput = fopen(filename , "rb" );
  if (finput==NULL) {fputs ("File error",stderr); exit (13);}
  
  fread(vector, sizeof(float), length, finput);
  for(i=0; i<length; i++){
    vector[i]=256*vector[i];
  }
  fclose(finput);
}

void read_image_pgm(unsigned char image[], char filename[], int imageWidth, int imageHeight)
{   /************************************************************************************
     * Function: void read_image_pgm(unsigned char image[], char filename[], int imageWidth, int imageHeight)
     * Input   : uchar array pointer for output result, char array with filename, int with with, int with height
     * Output  : uchar image array
     * Procedure: if image dimensions and layout pgm correct imare is read from file to image array
     ************************************************************************************/
  int grayMax;
  int PGM_HEADER_LINES=3;
  FILE* input;

  int headerLines = 1;
  int scannedLines= 0;
  long int counter =0;

  //read header strings
  char *lineBuffer = (char *) malloc(LINE_BUFFER_SIZE+1);
  char *split;
  char *format = (char *) malloc(LINE_BUFFER_SIZE+1);
  char P5[]="P5";
  char comments[LINE_BUFFER_SIZE+1];

  //open the input PGM file
  input=fopen(filename, "rb");

  //read the input PGM file header
  while(scannedLines < headerLines){
    fgets(lineBuffer, LINE_BUFFER_SIZE, input);
    //if not comments
    if(lineBuffer[0] != '#'){
      scannedLines += 1;
      //read the format
      if(scannedLines==1){
        split=strtok(lineBuffer, " \n");
        strcpy(format,split);
        if(strcmp(format,P5) == 0){
          //printf("FORMAT: %s\n",format);
          headerLines=PGM_HEADER_LINES;
        }
        else
        {
          printf("Only PGM P5 format is support. \n");
        }
      }
      //read width and height
      if (scannedLines==2)
      {
        split=strtok(lineBuffer, " \n");
        if(imageWidth == atoi(split)){ //check if width matches description
          //printf("WIDTH: %d, ", imageWidth);
        }
        else{
          printf("input frame has wrong width should be WIDTH: %d, ", imageWidth);
          exit(4);
        }
        split = strtok (NULL, " \n");
        if(imageHeight == atoi(split)){ //check if heigth matches description
          //printf("HEIGHT: %d\n", imageHeight);
        }
        else{
          printf("input frame has wrong height should be HEIGHT: %d, ", imageHeight);
          exit(4);
        }
      }
      // read maximum gray value
      if (scannedLines==3)
      {
        split=strtok(lineBuffer, " \n");
        grayMax = atoi(split);
        //printf("GRAYMAX: %d\n", grayMax);
      }
    }
    else
    {
      strcpy(comments,lineBuffer);
      //printf("comments: %s", comments);
    }
  }

  counter = fread(image, sizeof(unsigned char), imageWidth * imageHeight, input);
  //printf("pixels read: %d\n",counter);
        
  //close the input pgm file and free line buffer
  fclose(input);
  free(lineBuffer);
  free(format);
}

void write_image_pgm(unsigned char image[], const char filename[], int imageWidth, int imageHeight){
    FILE * output;
    output = fopen(filename, "wb");
    fprintf(output, "P5\n");
    fprintf(output, "%d %d\n255\n",imageWidth, imageHeight);
    fwrite(image, sizeof(unsigned char), imageWidth * imageHeight, output);
    fclose(output);
}
