#pragma once

#include <cmath>
#include <fstream>
#include <GL/freeglut.h>

using namespace std;


/**
 * An RGB triple
 */
class Color {
	
public:

	//------------
	// MEMBER DATA
	//------------

	/**
	 * The red channel
	 */
	GLdouble r;

	/**
	 * The green channel
	 */
	GLdouble g;
	
	/**
	 * The the blue channel
	 */
	GLdouble b;


	//-------------
	// CONSTRUCTORS
	//-------------

	/**
	 * Default constructor, initializes member data to 0s
	 */
	Color() : r( 0.0 ), g( 0.0 ), b( 0.0 ) {

		return;
	}


	/**
	 * Full constructor
	 * @param r The red channel
	 * @param g The green channel
	 * @param b The blue channel
	 */
	Color( GLdouble r, GLdouble g, GLdouble b ) : r( r ), g( g ), b( b ) {
		
		return;
	}

	/**
	 * Copy constructor
	 * @param other The color to copy
	 */
	Color( const Color & other ) : r( other.r ), g( other.g ), b( other.b ) {

		return;
	}

	
	//----------
	// OPERATORS
	//----------

	/**
	 * Assignment operator
	 * @param other The color to copy
	 * @return A copy of this color after assignment
	 */
	Color operator = ( const Color & other ) {

		r = other.r;
		g = other.g;
		b = other.b;

		return *this;
	}


	//----------
	// FUNCTIONS
	//----------

	/**
	 * Sets the member data to 0s
	 */
	void clear() {

		r = 0.0;
		g = 0.0;
		b = 0.0;
	}
};


/**
 * Assigns values to member data by reading from an input stream
 * @param stream The input file stream to read from
 * @param color The color to store the values in
 * @return The input stream (for cascading)
 */
ifstream & operator >> ( ifstream & stream, Color & color ) {

	stream >> color.r >> color.g >> color.b;

	return stream;
}
