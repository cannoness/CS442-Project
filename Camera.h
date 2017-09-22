#include <GL/freeglut.h>
#include "Vector.h"

/**
 * The number PI
 */
const GLdouble PI = 3.14159265358979;

/**
 * PI / 180 (for conversion between angle and radians)
 */
const GLdouble PI_DIV_180 = PI / 180.0;


/**
 * The camera for the scene
 */
class Camera {

public:

	//------------
	// MEMBER DATA
	//------------

	/**
	 * The location of the camera in 3D space
	 */
	Vector position;

	/**
	 * The direction that is "up" for teh camera
	 */
	Vector up;

	/**
	 * The direction that is "left" for the camera
	 */
	Vector left;

	/**
	 * The direction that is "forward" for the camera
	 */
	Vector forward;

	/**
	 * The distance from the position to the look at point (along the forward Vector)
	 */
	GLdouble focalLength;


	//-------------
	// CONSTRUCTORS
	//-------------

	/**
	 * Default constructor places camera at origin looking down -y axis, with +z up and +x left
	 */
	Camera() : position(), up( 0.0, 0.0, 1.0 ), left( 1.0, 0.0, 0.0 ), forward( 0.0, -1.0, 0.0 ), focalLength( 1.0 ) {

		return;
	}

	
	/**
	 * Full constructor
	 * @param position The location of the camera in 3D space
	 * @param up The direction that is "up" for the vector
	 * @param forward The direction that is "forward" for the camera
	 * @param focalLength The distance from the camera to the look-at point
	 */
	Camera( const Vector & position, const Vector & up, const Vector & forward, GLdouble focalLength ) : position( position ), up( up.unit() ), forward( forward.unit() ), left( ( forward ^ up ).unit() ), focalLength( focalLength ) {

		return;
	}


	/**
	 * Initializes the camera's member data using gluLookAt vectors
	 * @param lookFrom The location of the camera in 3D space
	 * @param lookAt The point the camera is focused on
	 * @param lookUp The direction that is "up" for the camera
	 */
	Camera( const Vector & lookFrom, const Vector & lookAt, const Vector & lookUp ) : position( lookFrom ), up( lookUp.unit() ), forward( ( lookAt - lookFrom ).unit() ), left( ( ( lookAt - lookFrom ) ^ lookUp ).unit() ), focalLength( ( lookAt - lookFrom ).length() ) {

		return;
	}


	//----------
	// OPERATORS
	//----------

	/**
	 * Assignment operator
	 * @param other The camera to copy
	 * @return A copy of the camera after assignment
	 */
	Camera operator = ( const Camera & other ) {

		position = other.position;
		up = other.up;
		forward = other.forward;
		left = other.left;
		focalLength = other.focalLength;

		return *this;
	}


	//----------
	// FUNCTIONS
	//----------

	/**
	 * Moves the camera along its forward vector
	 * @param amount The amount to move the vector
	 */
	void moveForward( GLdouble amount ) {

		position += forward * amount;
	}


	/**
	 * Moves the camera along its left vector
	 * @param amount The amount to move the vector
	 */
	void moveLeft( GLdouble amount ) {

		position += left * amount;
	}


	/**
	 * Moves the camera along its up vector
	 * @param amount The amount to move the vector
	 */
	void moveUp( GLdouble amount ) {

		position += up * amount;
	}


	/**
	 * Rotates the camera about its left axis
	 * @param angle The amount to rotate the camera
	 */
	void pitch( GLdouble angle ) {

		forward = ( forward * cos( angle * PI_DIV_180 ) + up * sin( angle * PI_DIV_180 ) ).unit();
		up = ( left ^ forward ).unit();
	}


	/**
	 * Rotates the camera about its up vector
	 * @param angle The amount to rotate the camera
	 */
	void yaw( GLdouble angle ) {

		forward = ( forward * cos( angle * PI_DIV_180 ) - left * sin( angle * PI_DIV_180 ) ).unit();
		left = ( forward ^ up ).unit();
	}


	/**
	 * Rotates the camera about its forward vector
	 * @param angle The amount to rotate the camera
	 */
	void roll( GLdouble angle ) {

		left = ( left * cos( angle * PI_DIV_180 ) + up * sin( angle * PI_DIV_180 ) ).unit();
		up = ( left ^ forward ).unit();
	}


	/**
	 * Moves the camera closer to or farther from its look-at point
	 * @param amount The amount to zoom
	 */
	void zoom( GLdouble amount ) {

		if ( focalLength - amount < 0.0 ) {

			return;
		}

		Vector lookAt = position + forward * focalLength;
		focalLength -= amount;
		position = lookAt - forward * focalLength;
	}


	/**
	 * Rotates the camera vertically about its look-at point
	 * @param angle The amount to rotate the camera
	 */
	void orbitUp( GLdouble angle ) {

		moveForward( focalLength );
		pitch( angle );
		moveForward( -focalLength );
	}


	/**
	 * Rotates the camera horizontally about its look-at point
	 * @param angle The amount to rotate the camera
	 */
	void orbitLeft( GLdouble angle ) {

		moveForward( focalLength );
		yaw( angle );
		moveForward( -focalLength );
	}


	/**
	 * Calls gluLookAt using the member data
	 */
	void updateGluLookAt() const {

		Vector lookAt = position + forward * focalLength;
		gluLookAt( position.x, position.y, position.z, lookAt.x, lookAt.y, lookAt.z, up.x, up.y, up.z );
	}
};
