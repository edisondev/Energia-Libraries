///
/// * Modified by Edison Developers on 16/02/2015 by adding the following functions:
///	*					void normalizeVec(float* vector)
///	*					float dotProduct(float* vector_one, float* vector_two)
///	*					void scalar_product(float* matrixIn, double scalar)
///	*			(Original downloaded from: http://forum.43oh.com/topic/3321-energia-library-matrix-for-launchpad-fraunchpad-and-stellarpad/)
/// * Created by Charlie Matlack on 12/18/10. Changed name.
/// * Original code by RobH45345 on Arduino Forums, taken from unknown source.
/// * Rei Vilo - Apr 15, 2011 - Added void MatrixMath::MatrixRowPrint(float* matrixA, uint8_t m, uint8_t n, String label){
/// * Rei Vilo - Jan 26, 2013 - Initial port to Energia
/// * Rick - Jan 27, 2013 - http://forum.43oh.com/topic/3275-is-there-an-energia-matrix-math-library-for-the-msp430-and-stellaris/?p=28604
/// * Rei Vilo - Feb 04, 2013 - Revised library with examples
///
///
#include "MatrixMath.h"

#define NR_END 1

MatrixMath::MatrixMath()
{
}

// Matrix Printing Routine
// Uses tabs to separate numbers under assumption printed float width won't cause problems
void MatrixMath::MatrixPrint(float* matrixA, uint8_t m, uint8_t n, String label){
  // A = input matrix (m x n)
  uint8_t i, j; 
  Serial.println(); 
  Serial.println(label); 
  for (i=0; i<m; i++){
    for (j=0; j<n; j++){
      Serial.print(matrixA[n*i+j]); 
      Serial.print("\t"); 
    }
    Serial.println(); 
  }
}

// avenue33 - Apr 15, 2011

void MatrixMath::MatrixRowPrint(float* matrixA, uint8_t m, uint8_t n, String label){
  // A = input matrix (m x n)
  uint8_t i, j; 
  Serial.print(label); 
  Serial.print("\t"); 
  for (j=0; j<n; j++){
    Serial.print(matrixA[n*m+j]); 
    Serial.print("\t"); 
  }
  Serial.println(); 
}


void MatrixMath::MatrixCopy(float* matrixA, uint8_t n, uint8_t m, float* matrixB)
{
  uint8_t i, j, k; 
  for (i=0; i<m; i++)
    for(j=0; j<n; j++)
    {
      matrixB[n*i+j] = matrixA[n*i+j]; 
    }
}

//Matrix Multiplication Routine
// C = A*B
void MatrixMath::MatrixMult(float* matrixA, float* matrixB, uint8_t m, uint8_t p, uint8_t n, float* matrixC)
{
  // A = input matrix (m x p)
  // B = input matrix (p x n)
  // m = number of rows in A
  // p = number of columns in A = number of rows in B
  // n = number of columns in B
  // C = output matrix = A*B (m x n)
  uint8_t i, j, k; 
  for (i=0; i<m; i++)
    for(j=0; j<n; j++)
    {
      matrixC[n*i+j]=0; 
      for (k=0; k<p; k++)
        matrixC[n*i+j]= matrixC[n*i+j]+matrixA[p*i+k]*matrixB[n*k+j]; 
    }
}


//Matrix Addition Routine
void MatrixMath::MatrixAdd(float* matrixA, float* matrixB, uint8_t m, uint8_t n, float* matrixC)
{
  // A = input matrix (m x n)
  // B = input matrix (m x n)
  // m = number of rows in A = number of rows in B
  // n = number of columns in A = number of columns in B
  // C = output matrix = A+B (m x n)
  uint8_t i, j; 
  for (i=0; i<m; i++)
    for(j=0; j<n; j++)
      matrixC[n*i+j]=matrixA[n*i+j]+matrixB[n*i+j]; 
}


//Matrix Subtraction Routine
void MatrixMath::MatrixSubtract(float* matrixA, float* matrixB, uint8_t m, uint8_t n, float* matrixC)
{
  // A = input matrix (m x n)
  // B = input matrix (m x n)
  // m = number of rows in A = number of rows in B
  // n = number of columns in A = number of columns in B
  // C = output matrix = A-B (m x n)
  uint8_t i, j; 
  for (i=0; i<m; i++)
    for(j=0; j<n; j++)
      matrixC[n*i+j] = matrixA[n*i+j] - matrixB[n*i+j]; 
}


//Matrix Transpose Routine
void MatrixMath::MatrixTranspose(float* matrixA, uint8_t m, uint8_t n, float* matrixC)

{
  // matrixA = input matrix (m x n)
  // m = number of rows in A
  // n = number of columns in A
  // matrixC = output matrix = the transpose of A (n x m)
  uint8_t i, j; 
  for (i=0; i<m; i++)
    for(j=0; j<n; j++)
      matrixC[m*j+i] = matrixA[n*i+j]; 
}


//Matrix Inversion Routine
// * This function inverts a matrix based on the Gauss Jordan method.
// * Specifically, it uses partial pivoting to improve numeric stability.
// * The algorithm is drawn from those presented in 
//	 NUMERICAL RECIPES: The Art of Scientific Computing.
// * The function returns 1 on success, 0 on failure.
// * NOTE: The argument is ALSO the result matrix, meaning the input matrix is REPLACED
uint8_t MatrixMath::MatrixInvert(float* matrixA, uint8_t n)
{
  // A = input matrix AND result matrix
  // n = number of rows = number of columns in A (n x n)
  uint8_t pivrow; 		// keeps track of current pivot row
  int8_t k, i, j; 		// k: overall index along diagonal; i: row index; j: col index
  uint8_t pivrows[n]; // keeps track of rows swaps to undo at end
  float tmp; 		// used for finding max value and making column swaps

  for (k = 0; k < n; k++)
  {
    // find pivot row, the row with biggest entry in current column
    tmp = 0; 
    for (i = k; i < n; i++) {
      if (abs(matrixA[i*n+k]) >= tmp)	{ // 'Avoid using other functions inside abs()?'
        tmp = abs(matrixA[i*n+k]); 
        pivrow = i; 
      }
    }

    // check for singular matrix
    if (matrixA[pivrow*n+k] == 0.0f)
    {
      Serial.println("Inversion failed due to singular matrix"); 
      return 0; 
    }

    // Execute pivot (row swap) if needed
    if (pivrow != k) {
      // swap row k with pivrow
      for (j = 0; j < n; j++) {
        tmp = matrixA[k*n+j]; 
        matrixA[k*n+j] = matrixA[pivrow*n+j]; 
        matrixA[pivrow*n+j] = tmp; 
      }
    }
    pivrows[k] = pivrow; 	// record row swap (even if no swap happened)

    tmp = 1.0f/matrixA[k*n+k]; 	// invert pivot element
    matrixA[k*n+k] = 1.0f; 		// This element of input matrix becomes result matrix

    // Perform row reduction (divide every element by pivot)
    for (j = 0; j < n; j++) {
      matrixA[k*n+j] = matrixA[k*n+j]*tmp; 
    }

    // Now eliminate all other entries in this column
    for (i = 0; i < n; i++) {
      if (i != k) {
        tmp = matrixA[i*n+k]; 
        matrixA[i*n+k] = 0.0f;  // The other place where in matrix becomes result mat
        for (j = 0; j < n; j++)
        {
          matrixA[i*n+j] = matrixA[i*n+j] - matrixA[k*n+j]*tmp; 
        }
      }
    }
  }

  // Done, now need to undo pivot row swaps by doing column swaps in reverse order
  for (k = n-1; k >= 0; k--) {
    if (pivrows[k] != k) {
      for (i = 0; i < n; i++) {
        tmp = matrixA[i*n+k];
        matrixA[i*n+k] = matrixA[i*n+pivrows[k]]; 
        matrixA[i*n+pivrows[k]] = tmp; 
      }
    }
  }
  return 1; 
}

float MatrixMath::DotProduct(float* vector_one, int length , float* vector_two){
  //This function computes a dotProduct of two equal sized vectors
  //Currently, there is no check implemented for equal-sizedness
  //Parameters:
  //  vector_one:  first input vector
  //  int size:    length of vector
  //  vector_two:  second input vector
  float out=0.0;
  for (uint8_t i=0; i<length;i++){
    out=out+vector_one[i]*vector_two[i];
  }
  return out;
}

void MatrixMath::NormalizeVector(float* MatIn, int n, int m, float* MatOut) {
  //This function normalizes the vector (but also works on matrices) so that norm|MatIn|=1
  //Input Parameters:
  //      MatIn:  Input Matrix/Vector
  //      n:      Number of rows
  //      m:      Number of Columns
  //      MatOut: Output Matrix
  float normVal=0.0;
  for (uint8_t i=0; i<m; i++){
    for (uint8_t j=0; j<n; j++){
      normVal=normVal+pow((double)MatIn[n*i+j],2.0);
    }
  }
  normVal=sqrt(normVal);
  Serial.println(normVal);
  if (normVal>1e-10){
    for (uint8_t i=0; i<m; i++){
      for (uint8_t j=0; j<n; j++){
        MatOut[n*i+j]=MatIn[n*i+j]/normVal;
      }
    }
  }else{  
	MatrixCopy((float*)MatIn, n, m, (float*)MatOut);   
  }
}


void MatrixMath::ScalarProduct(float* MatIn, int n, int m, double scalar,float* MatOut){
  //This function multiplies each element in the matrix by the scalar
  //Input Parameters:
  //      MatIn:  Input Matrix/Vector
  //      n:      Number of rows
  //      m:      Number of Columns
  //      scalar: Multiplication factor
  //      MatOut: Output Matrix  
  for (uint8_t i=0; i<m; i++){
    for (uint8_t j=0; j<n; j++){
      MatOut[n*i+j]=MatIn[n*i+j]*scalar;
    }
  }
  
}
