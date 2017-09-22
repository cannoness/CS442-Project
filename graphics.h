// Code to interface with the graphics from C.

/**
 * @param fileName The name of the planets file to read in order to initialize the planets list
 */
#ifdef __cplusplus
extern "C" {
#endif

void initParticleDisplay(int argc, char **argv);
void startParticles(void);
void enqueueParticle(double x, double y, double z, double r, double g, double b);
void finishParticles(void);
void updateStatistics(double deltaT);
void GraphicsDisplayStart();
void GraphicsDisplayFinish();
void GraphicsStart();
void GraphicsEnd();
void SimulationDisplay();
void SimulationDone();

#ifdef __cplusplus
}
#endif
