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

#define WIDTH 320*4
#define HEIGHT 240*4

using namespace std; 

vector<glm::vec3> vertices;
vector<ModelTriangle> faces;
vector<CanvasTriangle> canvasTriangles;
unordered_map<string, Colour> colourPalette;
vector<Colour> triangleColours;
float depthBuffer[WIDTH][HEIGHT];

struct camera {
    glm::vec3 position = glm::vec3(0,0,10);
    float focalLength = 10.0;
    glm::vec3 right = glm::vec3(1,0,0); glm::vec3 up = glm::vec3(0,1,0); glm::vec3 forward = glm::vec3(0,0,1);
    glm::mat3 orientation = glm::mat3(right, up, forward);
} camera;

glm::mat3 xRotationMatrix(float angle) {
    double pi = 2*acos(0.0);
    float angleRadian = angle * (pi / 180.0);
    glm::mat3 matrix = glm::mat3(
        1, 0, 0,
        0, cosf(angleRadian), sinf(angleRadian),
        0, -sinf(angleRadian), cosf(angleRadian)
    );
    return matrix;
}
glm::mat3 yRotationMatrix(float angle) {
    double pi = 2*acos(0.0);
    float angleRadian = angle * (pi / 180.0);
    glm::mat3 matrix = glm::mat3(
        cosf(angleRadian), 0, -sinf(angleRadian),
        0, 1, 0,
        sinf(angleRadian), 0, cosf(angleRadian)
    );
    return matrix;
}

void clearDB() {
    for(int i=0; i<WIDTH; i++){
		for(int j=0; j<HEIGHT; j++){
			depthBuffer[i][j] = 0;
		}
	}
}
glm::vec3 getVertexRelativeToCamera(glm::vec3 vertexPosition) {
	glm::vec3 cameraToVertex = camera.position - vertexPosition;
    glm::vec3 adjustedVector = cameraToVertex  * camera.orientation;
	return adjustedVector;
}
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength, DrawingWindow &window) {
	glm::vec3 relativePos = getVertexRelativeToCamera(vertexPosition);
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
void lineDrawingDepth(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour col){
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float depthDiff = to.depth - from.depth;
	float numberOfSteps = fmax(abs(xDiff), abs(yDiff))*2;
	float stepDepthDifference = depthDiff/numberOfSteps;
	vector<float> xValues = interpolateSingleFloats(from.x, to.x, numberOfSteps+1);
	vector<float> yValues = interpolateSingleFloats(from.y, to.y, numberOfSteps+1);
	uint32_t colour = (255<<24) + (int(col.red) << 16) +(int(col.green) << 8) +int(col.blue);
	for (float i=0; i<numberOfSteps; i++){
		int x = round(xValues[i]);
		int y = round(yValues[i]);
		float depth = from.depth+stepDepthDifference*i;
		if(depthBuffer[x][y] == 0){
			window.setPixelColour(x, y, colour);
			depthBuffer[x][y] = 1/depth;
		}
		else if(depthBuffer[x][y]>(1/depth)){
			window.setPixelColour(x, y, colour);
			depthBuffer[x][y] = 1/depth;
		}
	}
}
void stroked(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	lineDrawingDepth(window, triangle.v0(), triangle.v1(), colour);
	lineDrawingDepth(window, triangle.v0(), triangle.v2(), colour);
	lineDrawingDepth(window, triangle.v1(), triangle.v2(), colour);
}
void wireframeRender(DrawingWindow &window, glm::vec3 cameraPosition, float focalLength) {
	window.clearPixels();
	for (int i=0; i<faces.size(); i++) {
		Colour colour = Colour(255,255,255);
		vector<CanvasPoint> points;
		for (int j=0; j<faces[i].vertices.size(); j++) {
			CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, faces[i].vertices[j], focalLength, window);
			points.push_back(point);
		}
		CanvasTriangle face = CanvasTriangle(points[0], points[1], points[2]);
		points.clear();
		stroked(window, face, colour);
	}
}

CanvasTriangle sortTriangle(CanvasTriangle triangle) {
	if(triangle.v0().y < triangle.v1().y) {
		swap(triangle.v1(), triangle.v0());
	}
	if(triangle.v0().y < triangle.v2().y) {
		swap(triangle.v2(), triangle.v0());
	}
	if(triangle.v1().y < triangle.v2().y) {
		swap(triangle.v2(), triangle.v1());
	}
	return triangle;
}

CanvasPoint getMiddlePoint(CanvasTriangle triangle) {
	float proportion = (triangle.v1().y - triangle.v0().y)/(triangle.v2().y - triangle.v0().y);
	float proportionDepth = (triangle.v2().depth - triangle.v0().depth)/(triangle.v2().y - triangle.v0().y);
	float xCoord = triangle.v0().x + (proportion * (triangle.v2().x - triangle.v0().x));
	float yCoord = triangle.v1().y;
	float zCoord = triangle.v0().depth + (proportionDepth * (triangle.v1().y - triangle.v0().y));
	return CanvasPoint(xCoord, yCoord, zCoord);
}
vector<CanvasPoint> interpolateLine(CanvasPoint start, CanvasPoint end){
	vector<CanvasPoint> v;
	int numberOfValues = abs(start.y - end.y);
	float individualSeperation = (start.x - end.x)/numberOfValues;
	float individualDepthSepperation = -(start.depth - end.depth)/numberOfValues;
	if(start.y > end.y){
		for(int i=0; i<=numberOfValues; i++){ 
			v.push_back(CanvasPoint(round(start.x + i*(-individualSeperation)), round(start.y - i), (start.depth + i*(individualDepthSepperation))));
		}
	}
	else {
		for(int i=0; i<=numberOfValues; i++){ 
			v.push_back(CanvasPoint(round(start.x + i*(-individualSeperation)), round(start.y + i), (start.depth + i*(individualDepthSepperation))));
		}
	}
	return v;
}

void filled(DrawingWindow &window, CanvasTriangle tri, Colour col){
	tri = sortTriangle(tri);
	CanvasPoint newp = getMiddlePoint(tri);
	CanvasTriangle topTri = CanvasTriangle(tri.v0(), tri.v1(), newp);
	CanvasTriangle bottomTri = CanvasTriangle(tri.v2(), tri.v1(), newp);

	vector<CanvasPoint> side1Top = interpolateLine(topTri.v0(), topTri.v1());
	vector<CanvasPoint> side2Top = interpolateLine(topTri.v0(), topTri.v2());
	vector<CanvasPoint> side1Bottom = interpolateLine(bottomTri.v0(), bottomTri.v1());
	vector<CanvasPoint> side2Bottom = interpolateLine(bottomTri.v0(), bottomTri.v2());
	for(int i = 0; i< side1Top.size(); i++){
		lineDrawingDepth(window, side1Top[i], side2Top[i], col);
	}
	for(int i = 0; i< side1Bottom.size(); i++){
		lineDrawingDepth(window, side1Bottom[i], side2Bottom[i], col);
	}
	stroked(window, tri, col);
}

void rasterisedRender(DrawingWindow &window) {
    for(int i=0; i<canvasTriangles.size(); i++){
		filled(window, canvasTriangles[i], Colour(255,255,255));
	}
}

void lookAt() {
    glm::vec3 forward = camera.position;
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0),forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
    camera.orientation = glm::mat3(right, up, forward);
} 
void lookAt(glm::vec3 pointToLookAt) {
    glm::vec3 forward = camera.position - pointToLookAt;
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0),forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
    camera.orientation = glm::mat3(right, up, forward);
} 

void doStuff(DrawingWindow &window) {
    vector<CanvasPoint> canvasPoints;
	for(int i=0; i<faces.size(); i++){
		for(int j=0; j<3; j++){
			CanvasPoint temp = getCanvasIntersectionPoint(camera.position, faces[i].vertices[j], camera.focalLength, window);
			canvasPoints.push_back(temp);
		}
		triangleColours.push_back(faces[i].colour);
	}
	//puts the canvas points back into triangles but now with coords for the canvas
	for(int i=0; i<canvasPoints.size(); i=i+3){
		CanvasTriangle temp = CanvasTriangle(canvasPoints[i], canvasPoints[i+1], canvasPoints[i+2]);
		canvasTriangles.push_back(temp);
	}
}

void clearAll() {
    canvasTriangles.clear();
    triangleColours.clear();
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
		if (event.key.keysym.sym == SDLK_LEFT) {
            glm::vec3 left = glm::vec3(0.1, 0.0, 0.0);
            camera.position = camera.position + left;
            cout << "LEFT" << endl;
        } else if (event.key.keysym.sym == SDLK_RIGHT) {
            glm::vec3 right = glm::vec3(-0.1, 0.0, 0.0);
            camera.position = camera.position + right;
            cout << "RIGHT" << endl;
        } else if (event.key.keysym.sym == SDLK_UP) {
            glm::vec3 up = glm::vec3(0.0, -0.1, 0.0);
            camera.position = camera.position + up;
            cout << "UP" << endl;
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            glm::vec3 down = glm::vec3(0.0, 0.1, 0.0);
            camera.position = camera.position + down;
            cout << "DOWN" << endl;
        } else if (event.key.keysym.sym == SDLK_f) {
            glm::vec3 forward = glm::vec3(0.0, 0.0, -0.1);
            camera.position = camera.position + forward;
            cout << "FORWARD" << endl;
        } else if (event.key.keysym.sym == SDLK_b) {
            glm::vec3 backward = glm::vec3(0.0, 0.0, 0.1);
            camera.position = camera.position + backward;
            cout << "BACKWARD" << endl;
        } else if (event.key.keysym.sym == SDLK_x) {
            glm::mat3 rotationMatrix = xRotationMatrix(1.0);
            camera.position = camera.position * rotationMatrix;
            cout << "ROTATE X" << endl;
        } else if (event.key.keysym.sym == SDLK_y) {
            glm::mat3 rotationMatrix = yRotationMatrix(1.0);
            camera.position = camera.position * rotationMatrix;
            cout << "ROTATE Y" << endl;
        } else if (event.key.keysym.sym == SDLK_p) {
            glm::mat3 rotationMatrix = yRotationMatrix(1.0);
            camera.orientation = rotationMatrix * camera.orientation;
        } else if (event.key.keysym.sym == SDLK_t) {
            glm::mat3 rotationMatrix = xRotationMatrix(1.0);
            camera.orientation = rotationMatrix * camera.orientation;
        } else if (event.key.keysym.sym == SDLK_w) {
			clearAll();
			clearDB();
			window.clearPixels();
            wireframeRender(window, camera.position, camera.focalLength);
        } else if (event.key.keysym.sym == SDLK_r) {
			clearAll();
			clearDB();
            window.clearPixels();
			doStuff(window);
			rasterisedRender(window);
        } 
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
        // rasterisedRender(window, triangleColours);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}