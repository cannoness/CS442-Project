#include <GL/glut.h>
#include "Camera.h"
#include "graphics.h"

using namespace std;

//---------------
// CAMERA VECTORS
//---------------

Camera camera;
Camera originalCamera;

//---------------
// WINDOW DRAWING
//---------------

int calculations = 0;
GLint tmpParticles = 0;
GLint totalParticles = 0;
const GLint BORDER_WIDTH = 8;
const GLint TASKBAR_HEIGHT = 64;
GLint windowHeight = 0;
GLint windowWidth = 0;
GLboolean displayMenu = GL_FALSE;
GLboolean displayCPS = GL_TRUE;
char menu[] =
	"            CONTROLS\n"
	"------------------------------------\n"
	"         <ESC> - Quit\n"
	"          <F1> - Toggle this menu\n"
	"             f - Toggle Calculations Per Second\n"
	"             w - move camera forward\n"
	"             d - move camera backward\n"
	"             a - strafe camera left\n"
	"             d - strafe camera right\n"
	"             e - move camera up\n"
	"             c - move camera down\n"
	"             = - zoom in (move camera closer to look-at point\n"
	"             - - zoom out (move camera farther from look-at point\n"
	"  <hold shift> - 10x all keyboard camera movements\n"
	"         space - reorient view to original settings\n"
	"  <left mouse> - pan camera\n"
	" <right mouse> - orbit camera\n"
	"<middle mouse> - roll/zoom camera\n";


//--------
// PLANETS
//--------

GLdouble dt = 0.0;


//-------------
// KEY MOVEMENT
//-------------

GLdouble trackAmount = 1.0;


//---------------
// MOUSE MOVEMENT
//---------------

GLint prevX = -1;
GLint prevY = -1;
GLboolean active = GL_TRUE;
GLdouble panAmount = PI / 180.0;
const GLuint NO_MODE = 0;
const GLuint PAN_MODE = 1;
const GLuint ORBIT_MODE = 2;
const GLuint ROLL_MODE = 3;
GLuint mode = NO_MODE;



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

//-----------------
// HELPER FUNCTIONS
//-----------------

/**
 * Draws the coordinate axes in the scene (scaled by distance from look-from point to look-at point)
 */
void drawAxes() {

	glDisable( GL_LIGHTING );

	GLdouble magnitude = camera.focalLength / 3.0;

	glPushMatrix();
		glLineWidth( 1.0 );
		glBegin( GL_LINES );
			glColor3d( 1.0, 0.0, 0.0 );
			glVertex3d( -magnitude, 0.0, 0.0 );
			glVertex3d( magnitude, 0.0, 0.0 );

			glColor3d( 1.0, 1.0, 0.0 );
			glVertex3d( 0.0, -magnitude, 0.0 );
			glVertex3d( 0.0, magnitude, 0.0 );
			
			glColor3d( 0.0, 0.0, 1.0 );
			glVertex3d( 0.0, 0.0, -magnitude );
			glVertex3d( 0.0, 0.0, magnitude );
		glEnd();

		char label[] = "X-AXIS";
		glColor3d( 1.0, 0.0, 0.0 );
		drawText3D( magnitude, 0.0, 0.0, label );

		label[ 0 ] = 'Y';
		glColor3d( 1.0, 1.0, 0.0 );
		drawText3D( 0.0, magnitude, 0.0, label );

		label[ 0 ] = 'Z';
		glColor3d( 0.0, 0.0, 1.0 );
		drawText3D( 0.0, 0.0, magnitude, label );
	glPopMatrix();

	glEnable( GL_LIGHTING );
}


/**
 * If displayMenu is TRUE, then this function is called to display the CONTROLS menu
 */
void drawMenu() {
	
	glDisable( GL_LIGHTING );
	drawText2D( menu, 10, glutGet( GLUT_WINDOW_HEIGHT ) - 20 );
	glEnable( GL_LIGHTING );
}


/**
 * If displayCPS is TRUE, then this function is called to display the calculations-per-second AND the number of planets in the scene
 */
void drawCPS() {

	static clock_t then = clock();
	static GLdouble time = 0.0;
	
	char str[ 100 ];

	int calcs = calculations;

	if ( calcs >= 20 ) {

		clock_t now = clock();
		time = ( ( now - then ) / ( GLdouble ) ( CLOCKS_PER_SEC ) ) / ( GLdouble ) calcs;
		calculations = 0;
		then = now;
	}

	sprintf( str, "CPS: %2.5f\nParticles: %d\n", 1.0 / time, ( int ) totalParticles );

	glDisable( GL_LIGHTING );
	drawText2D( str, 10, 30 );
	glEnable( GL_LIGHTING );
}


//-------------------
// CALLBACK FUNCTIONS
//-------------------

void GraphicsDisplayStart() {
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        drawAxes();

        if ( displayMenu ) {

                drawMenu();
        }

        if ( displayCPS ) {

                drawCPS();
        }

}
void GraphicsDisplayFinish() {

        glutSwapBuffers();
        glFlush();

}
void GraphicsStart() {
        glutMainLoop();
}
void GraphicsEnd() {
	glutLeaveMainLoop();
}
void SimulationDisplay();
void SimulationDone();

/**
 * Called automatically when the window is resized
 * @param w The new width of the window
 * @param h The new height of the window
 */
void reshape( GLint w, GLint h ) {

	// Avoid dividing by 0
	if ( h == 0 ) {

		h = 1;
	}

	double ratio = w / ( double ) h;
	
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
    glViewport( 0, 0, w, h );
	gluPerspective( 45.0, ratio, 1.0, camera.focalLength * 2.0 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	camera.updateGluLookAt();
}


/**
 * Called automatically when a user presses a key
 * @param key The key that was pressed
 * @param x The x-coordinate of the mouse in window coordinates at the time of the key press
 * @param y The y-coordinate of the mosue in window coordinates at the time of the key press
 */
void keyboard( unsigned char key, GLint x, GLint y ) {

	GLdouble track = trackAmount;
	GLboolean moveLookAt = GL_TRUE;

	if ( glutGetModifiers() == GLUT_ACTIVE_SHIFT ) {

		track *= 10.0;
	}

	switch ( key ) {

		case 27: // ESC
			SimulationDone();
			break;

		case 'f': case 'F':

			displayCPS = !displayCPS;
			break;

		case 'w': case 'W':

			camera.moveForward( track );
			break;

		case 's': case 'S':

			camera.moveForward( -track );
			break;

		case 'a': case 'A':

			camera.moveLeft( -track );
			break;

		case 'd': case 'D':

			camera.moveLeft( track );
			break;

		case 'e': case 'E':

			camera.moveUp( track );
			break;

		case 'c': case 'C':

			camera.moveUp( -track );
			break;

		case '+': case '=':

			camera.zoom( track );
			break;

		case '-': case '_':

			camera.zoom( -track );
			break;

		case ' ':

			camera = originalCamera;
			break;
	}

	glLoadIdentity();
	camera.updateGluLookAt();
}


/**
 * Called automatically when a user presses a special key (like the function keys, for example)
 * @param key The special key that was pressed
 * @param x The x-coordinate of the mouse in window coordinates at the time of the key press
 * @param y The y-coordinate of the mouse in window coordinates at the time of the key press
 */
void specialFunc( GLint key, GLint x, GLint y ) {

	switch ( key ) {

		case GLUT_KEY_F1:

			displayMenu = !displayMenu;
			break;
	}
}


/**
 * Called automatically when the mouse is moved across the screen AND a button is currently pressed
 * @param x The new x-coordinate of the mouse in window coordinates
 * @param y The new y-coordinate of the mouse in window coordinates
 */
void motion( GLint x, GLint y ) {

	// If the mouse is outside the window, ignore the movement
	if ( !active ) {

		return;
	}

	if ( prevX == -1 || prevY == -1 ) {

		prevX = x;
		prevY = y;
	}

	GLdouble dx = ( GLdouble ) ( x - prevX );
	GLdouble dy = ( GLdouble ) ( y - prevY );

	// Left mouse button is pressed
	if ( mode == PAN_MODE ) {
		
		camera.yaw( -dx );
		camera.pitch( -dy );
	}

	// Right mouse button is pressed
	else if ( mode == ORBIT_MODE ) {
		
		camera.orbitLeft( dx );
		camera.orbitUp( dy );
	}

	// Middle mouse button is pressed
	else if ( mode == ROLL_MODE ) {

		camera.roll( dx );
		camera.zoom( dy * trackAmount );
	}

	glLoadIdentity();
	camera.updateGluLookAt();
	prevX = x;
	prevY = y;
}


/**
 * Called when the mouse moves across the screen AND no buttons are currently pressed
 * @param x The new x-coordinate of the mouse in window coordinates
 * @param y The new y-coordinate of the mouse in window coordinates
 */
void passiveMotion( GLint x, GLint y ) {

	prevX = x;
	prevY = y;
}


/**
 * Called automatically when the mouse enters or leaves the window
 * @param state The state of the mouse (entered or left the window)
 */
void entry( GLint state ) {

	if ( state == GLUT_ENTERED ) {
		
		active = GL_TRUE;
	}
	else {

		active = GL_FALSE;
	}
}


/**
 * Called automatically when a mouse button is pressed or released
 * @param button The button that was pressed
 * @param state The state of the mouse button (up or down)
 * @param x The x-coordinate of the mouse in window coordinates
 * @param y The y-coordinate of the mouse in window coordinates
 */
void mouse( GLint button, GLint state, GLint x, GLint y ) {

	if ( state == GLUT_UP ) {

		mode = 0;
	}
	else {

		if ( button == GLUT_LEFT_BUTTON ) {

			mode = PAN_MODE;
		}
		else if ( button == GLUT_RIGHT_BUTTON ) {

			mode = ORBIT_MODE;
		}
		else if ( button == GLUT_MIDDLE_BUTTON ) {

			mode = ROLL_MODE;
		}
	}
}


// Code to interface with the graphics from C.

/**
 * @param fileName The name of the planets file to read in order to initialize the planets list
 */
void initParticleDisplay(int argc, char **argv ) {
	totalParticles = 0;
	tmpParticles = 0;

        Vector lookFrom(0, 0, 100.0);
        Vector lookAt(0.0, 0.0, 0.0);
        Vector lookUp(0.0, 1.0, 0.0);

        glutInit( &argc, argv );
        glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA );
        glutInitWindowPosition( 0, 0 );

        windowWidth = glutGet( GLUT_SCREEN_WIDTH ) - BORDER_WIDTH;
        windowHeight = glutGet( GLUT_SCREEN_HEIGHT ) - TASKBAR_HEIGHT;
        glutInitWindowSize( windowWidth, windowHeight );
        glutCreateWindow( "Particles" );

        glutDisplayFunc( SimulationDisplay );
        glutReshapeFunc( reshape );
        glutIdleFunc( SimulationDisplay );
        glutKeyboardFunc( keyboard );
        glutSpecialFunc( specialFunc );
        glutMouseFunc( mouse );
        glutMotionFunc( motion );
        glutPassiveMotionFunc( passiveMotion );
        glutEntryFunc( entry );
        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	glEnable( GL_COLOR_MATERIAL );
	glColorMaterial( GL_FRONT, GL_AMBIENT_AND_DIFFUSE );

	// Somewhere in the initialization part of your program…
	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );
	 
	// Create light components
	GLfloat ambientLight[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
	GLfloat specularLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat positionLight[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	 
	// Assign created components to GL_LIGHT0
	glLightfv( GL_LIGHT0, GL_AMBIENT, ambientLight );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuseLight );
	glLightfv( GL_LIGHT0, GL_SPECULAR, specularLight );
	glLightfv( GL_LIGHT0, GL_POSITION, positionLight );

	glEnable( GL_DEPTH_TEST );
	glDepthMask( GL_TRUE );
	glEnable( GL_CULL_FACE );

	camera = Camera( lookFrom, lookAt, lookUp );
	originalCamera = camera;
}

void startParticles(void) {
	tmpParticles = 0;
	glPointSize( 2.0f );
	glBegin( GL_POINTS );
}

void enqueueParticle(double x, double y, double z,
	 	     double r, double g, double b) {
	glColor3d( r, g, b );
	float mcolor[] = { r, g, b, 1.0f };
	glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mcolor );
	glMaterialfv( GL_FRONT, GL_SPECULAR, mcolor );
	glVertex3d( x, y, z );
	tmpParticles++;
}

void finishParticles(void) {
	totalParticles = tmpParticles;
	calculations++;
	tmpParticles = 0;
	glEnd();
}

void updateStatistics(double deltaT) {
    dt = deltaT;
}

