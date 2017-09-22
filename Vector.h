#pragma once

#include <cmath>
#include <fstream>
#include <GL/freeglut.h>

using namespace std;


/**
 * A classic vector in 3D space
 */
class Vector {
	
public:

	//------------
	// MEMBER DATA
	//------------

	/**
	 * The x-coordinate of the vector
	 */
	GLdouble x;

	/**
	 * The y-coordinate of the vector
	 */
	GLdouble y;
	
	/**
	 * The z-coordinate of the vector
	 */
	GLdouble z;


	//-------------
	// CONSTRUCTORS
	//-------------

	/**
	 * Default constructor, initializes member data to 0s
	 */
	Vector() : x( 0.0 ), y( 0.0 ), z( 0.0 ) {

		return;
	}


	/**
	 * Full constructor
	 * @param x The x-coordinate of the vector
	 * @param y The y-coordiante of the vector
	 * @param z The z-coordiante of the vector
	 */
	Vector( GLdouble x, GLdouble y, GLdouble z ) : x( x ), y( y ), z( z ) {
		
		return;
	}

	/**
	 * Copy constructor
	 * @param other The vector to copy
	 */
	Vector( const Vector & other ) : x( other.x ), y( other.y ), z( other.z ) {

		return;
	}

	
	//----------
	// OPERATORS
	//----------

	/**
	 * Assignment operator
	 * @param other The vector to copy
	 * @return A copy of this vector after assignment
	 */
	Vector operator = ( const Vector & other ) {

		x = other.x;
		y = other.y;
		z = other.z;

		return *this;
	}


	/**
	 * Adds two vectors together
	 * @param other The vector to add to this one
	 * @return The sum of two vectors
	 */
	Vector operator + ( const Vector & other ) const {

		return Vector( x + other.x, y + other.y, z + other.z );
	}


	/**
	 * Adds two vectors together and stores the result in this vector
	 * @param other The vector to add to this one
	 * @return The sum of two vectors
	 */
	Vector operator += ( const Vector & other ) {

		x += other.x;
		y += other.y;
		z += other.z;

		return *this;
	}


	/**
	 * Subtracts the given vector from this one
	 * @param other The vector to subtract from this one
	 * @return The difference of the two vectors
	 */
	Vector operator - ( const Vector & other ) const {

		return Vector( x - other.x, y - other.y, z - other.z );
	}


	/**
	 * Subtracts the two vectors and stores the result in this vector
	 * @param other The vector to subtract from this one
	 * @return The difference of teh two vectors
	 */
	Vector operator -= ( const Vector & other ) {

		x -= other.x;
		y -= other.y;
		z -= other.z;

		return *this;
	}


	/**
	 * Creates the additive inverse of the vector
	 * @return The additive inverse of the vector
	 */
	Vector operator - () const {

		return Vector( -x, -y, -z );
	}

	/**
	 * Scales the vector by the given amount
	 * @param scalar The amount to scale the vector by
	 * @return The vector scaled by the given amount
	 */
	Vector operator * ( GLdouble scalar ) const {

		return Vector( x * scalar, y * scalar, z * scalar );
	}


	/**
	 * Computes the dot product of two vectors
	 * @param other The vector to dot with this one
	 * @return The dot product of the two vectors
	 */
	GLdouble operator * ( const Vector & other ) const {

		return x * other.x + y * other.y + z * other.z;
	}


	/**
	 * Scales the vector by the given amount and stores the result in this vector
	 * @param scalar The amount to scale the vector by
	 * @return The vector scaled by the given amount
	 */
	Vector operator *= ( GLdouble scalar ) {

		x *= scalar;
		y *= scalar;
		z *= scalar;

		return *this;
	}


	/**
	 * Shrinks the vector by the given amount
	 * @parma scalar The amount to shrink the vector by
	 * @return The vector shrunken by the given amount
	 */
	Vector operator / ( GLdouble scalar ) const {

		return Vector( x / scalar, y / scalar, z / scalar );
	}


	/**
	 * Shrinks the vector by the given amount and stores the result in this vector
	 * @param scalar The amount to shrink the vector by
	 * @return The vector shrunken by the given amount
	 */
	Vector operator /= ( GLdouble scalar ) {

		x /= scalar;
		y /= scalar;
		z /= scalar;

		return *this;
	}


	/**
	 * Computes the cross product of two vectors
	 * @param other The vector to cross with this one
	 * @return The cross product of two vectors
	 */
	Vector operator ^ ( const Vector & other ) const {

		return Vector( y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x );
	}


	/**
	 * Computes the cross product of two vectors and stores the result in this vector
	 * @param other The vector to cros with this one
	 * @return The cross product of two vectors
	 */
	Vector operator ^= ( const Vector & other ) {

		x = y * other.z - z * other.y;
		y = z * other.x - x * other.z;
		z = x * other.y - y * other.x;

		return *this;
	}


	//----------
	// FUNCTIONS
	//----------

	/**
	 * Computes the length of the vector
	 */
	GLdouble length() const {

		return sqrt( x * x + y * y + z * z );
	}


	/**
	 * Computes the square of the length of the vector
	 */
	GLdouble lengthSquared() const {

		return x * x + y * y + z * z;
	}


	/**
	 * Computes a unit vector in the same direction as this vector
	 * @return A unit-length vector in the same direction
	 */
	Vector unit() const {

		return Vector( x, y, z ) / length();
	}


	/**
	 * Normalizes this vector (length 1)
	 */
	void normalize() {

		GLdouble length = this->length();
		x /= length;
		y /= length;
		z /= length;
	}


	/**
	 * Sets the member data to 0s
	 */
	void clear() {

		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
};


/**
 * Assigns values to member data by reading from an input stream
 * @param stream The input file stream to read from
 * @param vector The vector to store the values in
 * @return The input stream (for cascading)
 */
ifstream & operator >> ( ifstream & stream, Vector & vector ) {

	stream >> vector.x >> vector.y >> vector.z;

	return stream;
}
