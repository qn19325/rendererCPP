#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <unordered_map>
#include <iostream>

#define WIDTH 320*2
#define HEIGHT 240*2

using namespace std;

vector<glm::vec3> vertices;
vector<ModelTriangle> faces;
unordered_map<string, Colour> colourPalette;

struct camera { 
	glm::vec3 position = glm::vec3(0.0, 0.0, 4.0);
    float focalLength = 2.0;
} camera;

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength, DrawingWindow &window) {
	glm::vec3 relativePos = vertexPosition - cameraPosition;
	float x = relativePos.x;
	float y = relativePos.y;
	float depth = relativePos.z;
	float u = 800 * (focalLength * (x / depth)) + (window.width / 2);
	float v = 800 * (focalLength * (y / depth)) + (window.height / 2);
	u = WIDTH - u;
	CanvasPoint point(u, v, depth);
	return point;
}
void pointCloud(DrawingWindow &window, glm::vec3 cameraPosition, float focalLength) {
	uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	for (int i=0; i<vertices.size(); i++)  {
		CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, vertices[i], focalLength, window);
		window.setPixelColour(point.x, point.y, colour);
	}
}
vector<float> interpolateSingleFloats(float from, float to, float numberOfValues) {
	vector<float> values;
	float difference = to - from;
	float stepSize = difference / (numberOfValues-1);
	for (int i=0; i<numberOfValues-1; i++) {
		values.push_back(from + stepSize*i);
	}
	return values;
}
void lineDrawing(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = fmax(abs(xDiff), abs(yDiff));
	vector<float> xValues = interpolateSingleFloats(from.x, to.x, numberOfSteps+1);
	vector<float> yValues = interpolateSingleFloats(from.y, to.y, numberOfSteps+1);
	uint32_t col = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + colour.blue;
	for (float i=0; i<numberOfSteps; i++) {
		int x = round(xValues[i]);
		int y = round(yValues[i]);
		window.setPixelColour(x, y, col);
	}
}
void stroked(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	lineDrawing(window, triangle.v0(), triangle.v1(), colour);
	lineDrawing(window, triangle.v0(), triangle.v2(), colour);
	lineDrawing(window, triangle.v1(), triangle.v2(), colour);
}
void wireframeRender(DrawingWindow &window, glm::vec3 cameraPosition, float focalLength) {
	for (int i=0; i<faces.size(); i++) {
		vector<CanvasPoint> points;
		for (int j=0; j<faces[i].vertices.size(); j++) {
			CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, faces[i].vertices[j], focalLength, window);
			points.push_back(point);
		}
		CanvasTriangle face = CanvasTriangle(points[0], points[1], points[2]);
		points.clear();
		stroked(window, face, Colour(255,255,255));
	}
}

void draw(DrawingWindow &window) {
	window.clearPixels();
}

void parseOBJ(string filePath) {
	float scaler = 0.001;
	ifstream fileIn;
	string line;
	string currentColour;
	fileIn.open(filePath);
	if (fileIn.is_open()) {
		while (getline(fileIn, line)) {
			if (line[0] == 'v' && line[1] == ' ') {
				vector<string> lineParts = split(line, ' ');
				glm::vec3 vertex = glm::vec3(stof(lineParts[1])*scaler,stof(lineParts[2])*scaler,stof(lineParts[3])*scaler);
				vertices.push_back(vertex);
			} else if (line[0] == 'f') {
				vector<string> lineParts = split(line, ' ');
				Colour colour = colourPalette[currentColour];
				vector<string> part1 = split(lineParts[1], '/');
				vector<string> part2 = split(lineParts[2], '/');
				vector<string> part3 = split(lineParts[3], '/');
				ModelTriangle face = ModelTriangle(vertices[stoi(part1[0])-1],vertices[stoi(part2[0])-1],vertices[stoi(part3[0])-1],colour);
				faces.push_back(face);
			} else if (line[0] == 'u') {
				vector<string> lineParts = split(line, ' ');
				currentColour = lineParts[1];
			}
		}
		fileIn.close();
	}
}

void parseMTL(string filePath) {
	ifstream fileIn;
	string line;
	Colour currentColour;
	string currentColourName;
	fileIn.open(filePath);
	if (fileIn.is_open()) {
		while (getline(fileIn, line)) {
			if(line[0] == 'n') {
				vector<string> lineParts = split(line, ' ');
				currentColourName = lineParts[1];
			} else if(line[0] == 'K') {
				vector<string> lineParts = split(line, ' ');
				int red = round(255 * stof(lineParts[1]));
				int green = round(255 * stof(lineParts[2]));
				int blue = round(255 * stof(lineParts[3]));
				currentColour = Colour(red,green,blue);
				colourPalette[currentColourName] = currentColour;
			}
		}
		fileIn.close();
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == 'w') window.clearPixels(); wireframeRender(window, camera.position, camera.focalLength);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	parseOBJ("build/logo/logo.obj");
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
