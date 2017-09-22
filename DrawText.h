#pragma once

#include <GL/freeglut.h>

/**
 * Draws the given text to the screen at the given screen coordinates (NOTE: in OpenGL, (0,0) is the bottom-left corner)
 * @param str The string to draw
 * @param x The x-coordinate of the left edge of the text in window coordinates
 * @param y The y-coordinate of the bottom of the text in window coordinates
 * @param r The red channel for the color of the text (default is 1.0)
 * @param g The green channel for the color of the text (default is 1.0)
 * @param b The blue channel for the color of the text (default is 1.0)
 */
void drawText2D( char * str, int x, int y, double r = 1.0, double g = 1.0, double b = 1.0 ) {

		glMatrixMode( GL_PROJECTION );
		glPushMatrix();
			glLoadIdentity();
			glOrtho( 0.0, glutGet( GLUT_WINDOW_WIDTH ), 0.0, glutGet( GLUT_WINDOW_HEIGHT ), -1.0, 1.0 );
			glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
					glLoadIdentity();
					glColor3d( r, g, b );
					glRasterPos2i( x, y );
					glutBitmapString( GLUT_BITMAP_8_BY_13, ( unsigned char * ) str );
				glPopMatrix();
			glMatrixMode( GL_PROJECTION );
		glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
}


/**
 * Draws the given text in world coordinates
 * @param x The x-coordinate for the text
 * @param y The y-coordinate for the text
 * @param z The z-coordinate for the text
 * @param str The string to draw
 */
void drawText3D( int x, int y, int z, char * str ) {

	glPushMatrix();
		glRasterPos3d( x, y, z );

		for ( char * c = str; *c != '\0'; ++c ) {

			glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, *c );
		}
	glPopMatrix();
}
