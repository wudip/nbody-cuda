#include <iostream>
#include <GLFW/glfw3.h>

#define WINDOW_WIDTH 640
#define WINDOW_HEIGHT 480
#define POINT_SIZE 20

using namespace std;

GLFWwindow* initWindow();
void drawRectangle(int x, int y);

int main(void) {
    cout << "Potato" << endl;
    GLFWwindow* window = initWindow();
    if (window == NULL) {
        cerr << "Cannot initialize GLFW library for displaying the output in a window." << endl;
        return -1;
    }

    glfwMakeContextCurrent(window);

    /* Some stuff to make the rectangle appear */
    glViewport(0.0f, 0.0f, WINDOW_WIDTH, WINDOW_HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, 0, 1);
    glMatrixMode(GL_MODELVIEW);

    int x = 50;
    int y = 30;
    int xDirection = 1;
    int yDirection = -1;

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window)) {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        drawRectangle(x, y);
        
        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();

        /*
         * Update coords
         */
        if (x <= POINT_SIZE / 2 || x >= WINDOW_WIDTH - 1 - POINT_SIZE / 2) {
          xDirection *= -1;
        }
        if (y <= POINT_SIZE / 2 || y >= WINDOW_HEIGHT - 1 - POINT_SIZE / 2) {
          yDirection *= -1;
        }
        x += xDirection;
        y += yDirection;
    }

    glfwTerminate();
    return 0;

    cout << "I like trains!" << endl;
    return 0;
}

/**
  * Creates GLFW window
  * @return GLFW window or NULL in case operation was not successful
  */
GLFWwindow* initWindow() {
  return glfwInit() ? glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "I like trains!", NULL, NULL) : NULL;
}

void drawRectangle(int x, int y) {
    GLfloat pointVertex[] = { (float)x, (float)y };
    glEnableClientState(GL_VERTEX_ARRAY);
    glPointSize(POINT_SIZE);
    glVertexPointer(2, GL_FLOAT, 0, pointVertex);
    glDrawArrays(GL_POINTS, 0, 1);
    glDisableClientState(GL_VERTEX_ARRAY);
}