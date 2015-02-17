///
/// @file	MatrixMath.h
/// @brief	Library header
/// @details	Matrix library
/// @n	
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
/// * Edison Developers - Feb 16, 2015 - Addition of the following functions:
///	*					void normalizeVec(float* vector)
///	*					float dotProduct(float* vector_one, float* vector_two)
///	*					void scalar_product(float* matrixIn, double scalar
///
/// @n
/// @n @a	Developed with [embedXcode](http://embedXcode.weebly.com)
/// @n
/// @author	Rei VILO
/// @author	embedXcode.weebly.com
/// @date	Feb 03, 2013
/// @version	103
///
/// @copyright	© Rei VILO, 2013
/// @copyright	CC = BY NC SA
///
/// @see	ReadMe.txt for references
#ifndef MatrixMath_h
#define MatrixMath_h
#include "Energia.h"
class MatrixMath
{
public:
	
	
	
	MatrixMath();
		
	void MatrixPrint(float* matrixA, uint8_t m, uint8_t n, String label);
	// vilo.rei - Apr 15, 2011	
	
	void MatrixRowPrint(float* matrixA, uint8_t m, uint8_t n, String label);
		
	void MatrixCopy(float* matrixA, uint8_t n, uint8_t m, float* B);
			
	void MatrixMult(float* matrixA, float* matrixB, uint8_t m, uint8_t p, uint8_t n, float* matrixC);
		
	void MatrixAdd(float* matrixA, float* matrixB, uint8_t m, uint8_t n, float* matrixC);

	void MatrixSubtract(float* matrixA, float* matrixB, uint8_t m, uint8_t n, float* matrixC);
	
	void MatrixTranspose(float* matrixA, uint8_t m, uint8_t n, float* matrixC);
	
	uint8_t MatrixInvert(float* matrixA, uint8_t n);
	
	float DotProduct(float* vector_one, int length , float* vector_two);
	//EdisonDev - Feb 16, 2015
	
	void NormalizeVector(float* MatIn, int n, int m, float* MatOut);
	//EdisonDev - Feb 16, 2015
	
	void ScalarProduct(float* MatIn, int n, int m, double scalar,float* MatOut);
	//EdisonDev - Feb 16, 2015
};
#endif
