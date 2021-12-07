#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <ModelTriangle.h>
#include <unordered_map>
#include <RayTriangleIntersection.h>
#include <math.h>

#define WIDTH 320*3
#define HEIGHT 240*3

using namespace std;

string model;
vector<glm::vec3> vertices;
vector<glm::vec3> normals;
vector<ModelTriangle> faces;
vector<CanvasTriangle> canvasTriangles;
unordered_map<string, Colour> colourPalette;
vector<Colour> triangleColours; 
float depthBuffer[WIDTH][HEIGHT];
string renderMode;
glm::vec3 lightPosition;
struct camera {
    glm::vec3 position;
    float focalLength;
} camera;

void clearDepthBuffer() {
    for(int i=0; i<WIDTH; i++){
		for(int j=0; j<HEIGHT; j++){
			depthBuffer[i][j] = 0;
		}
	}
}

glm::vec3 getVertexRelativeToCamera(glm::vec3 vertexPosition) {
	glm::vec3 relative = vertexPosition - camera.position;
	return relative;
}
CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, DrawingWindow &window) {
	glm::vec3 relativePos = getVertexRelativeToCamera(vertexPosition);
	float x = relativePos.x;
	float y = relativePos.y;
	float depth = relativePos.z;
	float u = 400 * (camera.focalLength * (x / depth)) + (window.width / 2);
	float v = 400 * (camera.focalLength * (y / depth)) + (window.height / 2);
	u = WIDTH - u;
	CanvasPoint point(u, v, depth);
	return point;
}
void pointCloud(DrawingWindow &window) {
    window.clearPixels();
	uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	for (int i=0; i<vertices.size(); i++)  {
		CanvasPoint point = getCanvasIntersectionPoint(vertices[i], window);
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
void strokedNoDepth(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	lineDrawing(window, triangle.v0(), triangle.v1(), colour);
	lineDrawing(window, triangle.v0(), triangle.v2(), colour);
	lineDrawing(window, triangle.v1(), triangle.v2(), colour);
}
void strokedDepth(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	lineDrawingDepth(window, triangle.v0(), triangle.v1(), colour);
	lineDrawingDepth(window, triangle.v0(), triangle.v2(), colour);
	lineDrawingDepth(window, triangle.v1(), triangle.v2(), colour);
}

void wireFrame(DrawingWindow &window) {
    window.clearPixels();
	for (int i=0; i<faces.size(); i++) {
		vector<CanvasPoint> points;
		for (int j=0; j<3; j++) {
			CanvasPoint point = getCanvasIntersectionPoint(faces[i].vertices[j], window);
			points.push_back(point);
		}
		CanvasTriangle face = CanvasTriangle(points[0], points[1], points[2]);
		points.clear();
		strokedNoDepth(window, face, Colour(255,255,255));
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
vector<CanvasPoint> interpolateLine(CanvasPoint start, CanvasPoint end, int side){
	vector<CanvasPoint> v;
	int numberOfValues = abs(start.y - end.y);
	float individualSeperation = (start.x - end.x)/numberOfValues;
	float individualDepthSepperation = -(start.depth - end.depth)/numberOfValues;
	if(side == 0){
		for(int i=0; i<=numberOfValues; i++){ 
			v.push_back(CanvasPoint(roundf(start.x + i*(-individualSeperation)), start.y - i, (start.depth + i*(individualDepthSepperation))));
		}
	}
	else if(side == 1){
		for(int i=0; i<=numberOfValues; i++){ 
			v.push_back(CanvasPoint(roundf(start.x + i*(-individualSeperation)), start.y + i, (start.depth + i*(individualDepthSepperation))));
		}
	}
	return v;
}

void interpolateSidesAndFill(CanvasTriangle tri, DrawingWindow &window, Colour col, int side){
	vector<CanvasPoint> side1 = interpolateLine(tri.v0(), tri.v1(), side);
	vector<CanvasPoint> side2 = interpolateLine(tri.v0(), tri.v2(), side);
	for(int i = 0; i< side1.size(); i++){
		lineDrawingDepth(window, side1[i], side2[i], col);
	}
	strokedDepth(window, tri, col);
}
void rasterise(DrawingWindow &window){
    window.clearPixels();
    clearDepthBuffer();
    canvasTriangles.clear();
    if (canvasTriangles.size() == 0) {
        vector<CanvasPoint> canvasPoints;
        for(int i=0; i<faces.size(); i++) {
            for(int j=0; j<3; j++){
                CanvasPoint temp = getCanvasIntersectionPoint(faces[i].vertices[j], window);
                canvasPoints.push_back(temp);
            }
            triangleColours.push_back(faces[i].colour);
        }
        for(int i=0; i<canvasPoints.size(); i=i+3){
            CanvasTriangle temp = CanvasTriangle(canvasPoints[i], canvasPoints[i+1], canvasPoints[i+2]);
            canvasTriangles.push_back(temp);
        }
    }
    for (int i=0; i<canvasTriangles.size(); i++) {
        Colour col = triangleColours[i];
        CanvasTriangle tri = canvasTriangles[i];
        tri = sortTriangle(tri);
        CanvasPoint newp = getMiddlePoint(tri);
        CanvasTriangle topTri = CanvasTriangle(tri.v0(), tri.v1(), newp);
        CanvasTriangle bottomTri = CanvasTriangle(tri.v2(), tri.v1(), newp);
        interpolateSidesAndFill(topTri, window, col, 0);
        interpolateSidesAndFill(bottomTri, window, col, 1);
    }

}

RayTriangleIntersection getClosestIntersectionRayTrace(glm::vec3 rayDirection, glm::vec3 rayStart, int index) {
	float closestSoFar = FLT_MAX;
	RayTriangleIntersection intersection;
	for (size_t i = 0; i < faces.size(); i++) {
		ModelTriangle triangle = faces[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = rayStart - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution[0]; float u = possibleSolution[1]; float v = possibleSolution[2];
		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((u + v) <= 1.0) && (t > 0 && t < closestSoFar) && i != index) {
			closestSoFar = t;
			intersection.intersectionPoint = rayStart + (rayDirection * t);
			intersection.distanceFromCamera = t;
			intersection.intersectedTriangle = triangle;
			intersection.triangleIndex = i;
		}
	}
	return intersection;
}

bool rayCanReachLight(glm::vec3 point, int index) {
	glm::vec3 rayDirection = lightPosition - point;
	RayTriangleIntersection intersection = getClosestIntersectionRayTrace(rayDirection, point, index);
	if (intersection.triangleIndex != -1) {
		float distanceToIntersection = glm::l2Norm(intersection.intersectionPoint, point);
		float distanceToLight = glm::l2Norm(lightPosition, point);
		if (distanceToIntersection > distanceToLight) return true;
		else return false;
	}
	return true;
}

void raytrace(DrawingWindow& window) {
	window.clearPixels();
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			float xRayDirection = (x - WIDTH / 2.0) ;
			float yRayDirection = y - HEIGHT / 2.0;
			xRayDirection = xRayDirection / 400;
			yRayDirection = -yRayDirection / 400;
			glm::vec3 imagePlanePoint(xRayDirection, yRayDirection, camera.position.z - camera.focalLength);
			glm::vec3 rayDirection = imagePlanePoint - camera.position;
			RayTriangleIntersection intersection = getClosestIntersectionRayTrace(rayDirection, camera.position, -1);
			if (intersection.triangleIndex != -1) {
				Colour colour;
                glm::vec3 lightRay = glm::vec3(lightPosition - intersection.intersectionPoint);
				// proximity
                float proximity = 1 / (2 * M_PI * glm::pow(glm::length(lightRay), 2));
				if (proximity > 1) proximity = 1;
                // angle of incidence
				float angleOfIncidence = glm::dot(glm::normalize(intersection.intersectedTriangle.faceNormal), glm::normalize(lightRay));
				if (angleOfIncidence < 0) angleOfIncidence = 0;
				// specular
				int n = 10;
				glm::vec3 lightToPoint = glm::normalize(intersection.intersectionPoint - lightPosition);
				glm::vec3 normalizedNormal = glm::normalize(intersection.intersectedTriangle.faceNormal);
				glm::vec3 rayReflection = glm::normalize(lightToPoint - normalizedNormal * 2.0f * glm::dot(lightToPoint, normalizedNormal)); 
				glm::vec3 toCamera = glm::normalize(camera.position - intersection.intersectionPoint);
				float specularLight = glm::dot(toCamera,rayReflection);
				if (specularLight < 0) specularLight = 0;
				specularLight = pow(specularLight, n);
				
				float brightness = (proximity + angleOfIncidence + specularLight) / 3;
				if ((brightness+0.2) < 1) brightness += 0.2;
				else brightness = 1.0;
				if (rayCanReachLight(intersection.intersectionPoint,intersection.triangleIndex)) {
					colour = intersection.intersectedTriangle.colour;
					uint32_t col = (255 << 24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + int(colour.blue*brightness);
					window.setPixelColour(x, y, col);
				}
				else {
					brightness = 0.2;
					colour = intersection.intersectedTriangle.colour;
					uint32_t col = (255 << 24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + int(colour.blue*brightness);
					window.setPixelColour(x, y, col);
				}
			}
		}
	}
}

RayTriangleIntersection getClosestIntersectionShading(glm::vec3 rayDirection, glm::vec3 rayStart) {
	float closestSoFar = FLT_MAX;
	int index = -1;
	RayTriangleIntersection intersection;
	intersection.triangleIndex = -1;
 	for(float i=0; i < faces.size(); i++){
		ModelTriangle triangle = faces[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = rayStart - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution[0]; float u = possibleSolution[1]; float v = possibleSolution[2];
		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((u + v) <= 1.0) && ((t > 0 && t < closestSoFar) || index == -1)) {
			closestSoFar = t; index = i;
			intersection.distanceFromCamera = t; 
			intersection.intersectedTriangle = triangle; 
			intersection.triangleIndex = i;
			intersection.intersectionPoint = triangle.vertices[0] + (e0*u) + (e1*v);	
			intersection.baryCoords[1] = u; intersection.baryCoords[2] = v; intersection.baryCoords[0] = 1 - (u + v);	
		}
	}
	return intersection;
}

void gouraud(DrawingWindow &window) {
	window.clearPixels();	
	for(float y = 0; y < HEIGHT; y++){
		for(float x = 0; x < WIDTH; x++){
			float xRayDirection = (x - WIDTH / 2.0) ;
			float yRayDirection = y - HEIGHT / 2.0;
			xRayDirection = xRayDirection / 400;
			yRayDirection = -yRayDirection / 400;
			glm::vec3 imagePlanePoint(xRayDirection, yRayDirection, camera.focalLength - camera.position.z );
			glm::vec3 rayDirection = imagePlanePoint - camera.position;
			RayTriangleIntersection intersection = getClosestIntersectionShading(rayDirection, camera.position);
			if(intersection.triangleIndex != -1) {
				ModelTriangle face = faces[intersection.triangleIndex];
				for(int i = 0; i < 3; i++) {
					glm::vec3 normal = face.vertexNormals[i];
					glm::vec3 lightRay =  glm::normalize(lightPosition - intersection.intersectionPoint);
        			float aoi = glm::dot(normal, lightRay);
        			if (aoi<0.2) aoi=0.2;
					// NEED SPECULAR
					int n = 256;
					glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - lightPosition);
					normal = glm::normalize(normal);
					glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2.0f*normal*(glm::dot(vectorOfLightToPoint, normal))));
					glm::vec3 vectorToCamera = glm::normalize((camera.position - intersection.intersectionPoint));
					float specBrightness = glm::dot(vectorToCamera, vectorOfReflection);
					if (specBrightness < 0) specBrightness = 0;
					specBrightness = pow(specBrightness, n);
					float brightness = (aoi + specBrightness)/2;
					face.vertexBrightness[i] = brightness;
				}
				float brightness = (face.vertexBrightness[0]*intersection.baryCoords[0]) + (face.vertexBrightness[1]*intersection.baryCoords[1]) + (face.vertexBrightness[2]*intersection.baryCoords[2]);
                Colour colour = Colour(255,0,0);
				uint32_t col = (255 << 24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + int(colour.blue*brightness);
				window.setPixelColour(x, y, col);
			}
		}
	}
}

array<float, 3> getVertexBrightnesses(RayTriangleIntersection intersection, glm::vec3 cameraPosition){
	array<float, 3> returnArray;
	for(int i = 0; i<3; i++){
		glm::vec3 vectorFromLight = glm::normalize(lightPosition - intersection.intersectionPoint);
		float aoiBrightness = glm::dot(intersection.intersectedTriangle.vertexNormals[i], vectorFromLight);
		
		int n = 100;
		glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - lightPosition);
		glm::vec3 normal = glm::normalize(intersection.intersectedTriangle.vertexNormals[i]);
		glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2.0f*normal*(glm::dot(vectorOfLightToPoint, normal))));
		glm::vec3 vectorToCamera = glm::normalize((cameraPosition - intersection.intersectionPoint));
		float specBrightness = glm::dot(vectorToCamera, vectorOfReflection);
		if (specBrightness < 0){
			specBrightness = 0;
		}
		specBrightness = pow(specBrightness, n);
		
		float vertexBrightness = (specBrightness + aoiBrightness)/2;
		if (vertexBrightness<0){
			vertexBrightness = 0;
		}
		if (vertexBrightness>1){
			vertexBrightness = 1;
		}
		returnArray[i] = vertexBrightness;
	}
	return returnArray;
}

void phong(DrawingWindow &window) {
	window.clearPixels();	
	for(float y = 0; y < HEIGHT; y++){
		for(float x = 0; x < WIDTH; x++){
			float xRayDirection = (x - WIDTH / 2.0) ;
			float yRayDirection = y - HEIGHT / 2.0;
			xRayDirection = xRayDirection / 400;
			yRayDirection = -yRayDirection / 400;
			glm::vec3 imagePlanePoint(xRayDirection, yRayDirection, camera.focalLength - camera.position.z );
			glm::vec3 rayDirection = imagePlanePoint - camera.position;
			RayTriangleIntersection intersection = getClosestIntersectionShading(rayDirection, camera.position);
			if(intersection.triangleIndex != -1) {
				glm::vec3 pixelNormal = glm::normalize(glm::normalize(intersection.intersectedTriangle.vertexNormals[0])*intersection.baryCoords[0] + glm::normalize(intersection.intersectedTriangle.vertexNormals[1])*intersection.baryCoords[1] + glm::normalize(intersection.intersectedTriangle.vertexNormals[2])*intersection.baryCoords[2]);

				glm::vec3 vectorFromLight = glm::normalize(lightPosition - intersection.intersectionPoint);
				float aoiBrightness = glm::dot(pixelNormal, vectorFromLight);
				if(aoiBrightness<0){
					aoiBrightness = 0;
				}

				int n = 256;
				glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - lightPosition);
				glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*pixelNormal*(glm::dot(vectorOfLightToPoint, pixelNormal))));
				glm::vec3 vectorToCamera = glm::normalize((camera.position - intersection.intersectionPoint));
				float specBrightness = glm::dot(vectorToCamera, vectorOfReflection);
				if (specBrightness < 0){
					specBrightness = 0;
				}
				specBrightness = pow(specBrightness, n);
				
				float pixelBrightness = (specBrightness + aoiBrightness);
				if (pixelBrightness<0){
					pixelBrightness = 0;
				}
				if (pixelBrightness>1){
					pixelBrightness = 1;
				}
				if (pixelBrightness + 0.1 <= 1.0) pixelBrightness += 0.1;
				else  pixelBrightness = 1.0;
				uint32_t col = (255<<24) + (int(255*pixelBrightness) << 16) + (int(0*pixelBrightness) << 8) + (int(0*pixelBrightness));
				window.setPixelColour(x, y, col);
			}

		}
	}
}

void reflection(DrawingWindow& window) {
	window.clearPixels();
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			float xRayDirection = (x - WIDTH / 2.0) ;
			float yRayDirection = y - HEIGHT / 2.0;
			xRayDirection = xRayDirection / 400;
			yRayDirection = -yRayDirection / 400;
			glm::vec3 imagePlanePoint(xRayDirection, yRayDirection, camera.position.z - camera.focalLength);
			glm::vec3 rayDirection = imagePlanePoint - camera.position;
			RayTriangleIntersection intersection = getClosestIntersectionRayTrace(rayDirection, camera.position, -1);
			if (intersection.intersectedTriangle.reflective == true) {
				glm::vec3 reflectionVector = rayDirection - (2.0f * intersection.intersectedTriangle.faceNormal * glm::dot(rayDirection, intersection.intersectedTriangle.faceNormal));
				RayTriangleIntersection reflectionIntersection = getClosestIntersectionRayTrace(glm::normalize(reflectionVector), intersection.intersectionPoint + 0.0001f * intersection.intersectedTriangle.faceNormal, -1);
				intersection = reflectionIntersection;
			}
			if (intersection.triangleIndex != -1) {
				Colour colour;
                glm::vec3 lightRay = glm::vec3(lightPosition - intersection.intersectionPoint);
				// proximity
                float proximity = 1 / (2 * M_PI * glm::pow(glm::length(lightRay), 2));
				if (proximity > 1) proximity = 1;
                // angle of incidence
				float angleOfIncidence = glm::dot(glm::normalize(intersection.intersectedTriangle.faceNormal), glm::normalize(lightRay));
				if (angleOfIncidence < 0) angleOfIncidence = 0;
				// specular
				int n = 10;
				glm::vec3 lightToPoint = glm::normalize(intersection.intersectionPoint - lightPosition);
				glm::vec3 normalizedNormal = glm::normalize(intersection.intersectedTriangle.faceNormal);
				glm::vec3 rayReflection = glm::normalize(lightToPoint - normalizedNormal * 2.0f * glm::dot(lightToPoint, normalizedNormal)); 
				glm::vec3 toCamera = glm::normalize(camera.position - intersection.intersectionPoint);
				float specularLight = glm::dot(toCamera,rayReflection);
				if (specularLight < 0) specularLight = 0;
				specularLight = pow(specularLight, n);
				
				float brightness = (proximity + angleOfIncidence + specularLight) / 3;
				if ((brightness+0.2) < 1) brightness += 0.2;
				else brightness = 1.0;
				if (rayCanReachLight(intersection.intersectionPoint,intersection.triangleIndex)) {
					colour = intersection.intersectedTriangle.colour;
					uint32_t col = (255 << 24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + int(colour.blue*brightness);
					window.setPixelColour(x, y, col);
				}
				else {
					brightness = 0.2;
					colour = intersection.intersectedTriangle.colour;
					uint32_t col = (255 << 24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + int(colour.blue*brightness);
					window.setPixelColour(x, y, col);
				}
			}
		}
	}
}

void clearAll() {
	vertices.clear();
	normals.clear();
	faces.clear();
	canvasTriangles.clear();
	colourPalette.clear();
	triangleColours.clear(); 
	clearDepthBuffer();
}

void parseOBJ(string filePath, float scaler) {
	ifstream fileIn;
	string line;
	string currentColour;
	bool reflective = false;
	fileIn.open(filePath);
	if (fileIn.is_open()) {
		while (getline(fileIn, line)) {
			if (line[0] == 'v') {
				vector<string> lineParts = split(line, ' ');
				if (line[1] == 'n') {
					glm::vec3 normal = glm::vec3(stof(lineParts[1]),stof(lineParts[2]),stof(lineParts[3]));
					normals.push_back(normal);
				} else {
					glm::vec3 vertex = glm::vec3(stof(lineParts[1])*scaler,stof(lineParts[2])*scaler,stof(lineParts[3])*scaler);
					vertices.push_back(vertex);
				}
			} else if (line[0] == 'f') {
				vector<string> lineParts = split(line, ' ');
				vector<string> splitLineParts1 = split(lineParts[1], '/');
				vector<string> splitLineParts2 = split(lineParts[2], '/');
				vector<string> splitLineParts3 = split(lineParts[3], '/');
				Colour colour = colourPalette[currentColour];
                glm::vec3 v0 = vertices[stoi(splitLineParts1[0])-1];
                glm::vec3 v1 = vertices[stoi(splitLineParts2[0])-1];
                glm::vec3 v2 = vertices[stoi(splitLineParts3[0])-1];
				ModelTriangle face = ModelTriangle(v0, v1, v2, colour);
                face.faceNormal = glm::normalize(glm::cross((v1-v0),(v2-v0)));
				if (normals.size() !=  0) {
					face.vertexNormals[0] = normals[stoi(splitLineParts1[2])-1];
					face.vertexNormals[1] = normals[stoi(splitLineParts2[2])-1];
					face.vertexNormals[2] = normals[stoi(splitLineParts3[2])-1];
				}
				face.reflective = reflective;
				faces.push_back(face);
			} else if (line[0] == 'u') {
				vector<string> lineParts = split(line, ' ');
				currentColour = lineParts[1];
			} else if (line[0] == 'r') {
                vector<string> lineParts = split(line, ' ');
                if (lineParts[1] == "f") reflective = false;
                else if (lineParts[1] == "t") reflective = true;
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

void draw(DrawingWindow &window) {
    if (renderMode == "pointCloud") pointCloud(window);
    else if (renderMode == "wireFrame") wireFrame(window);
    else if (renderMode == "rasterise") rasterise(window);
    else if (renderMode == "raytrace") raytrace(window);
	else if (renderMode == "gouraud") gouraud(window);
	else if (renderMode == "phong") phong(window);
	else if (renderMode == "reflection") reflection(window);
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) {
            glm::vec3 left = glm::vec3(0.1, 0.0, 0.0);
            camera.position = camera.position + left;
            cout << "LEFT" << endl; 
            window.renderFrame(); } 
        else if (event.key.keysym.sym == SDLK_RIGHT) {
            glm::vec3 right = glm::vec3(-0.1, 0.0, 0.0);
            camera.position = camera.position + right;
            cout << "RIGHT" << endl; } 
        else if (event.key.keysym.sym == SDLK_UP) {
            glm::vec3 up = glm::vec3(0.0, -0.1, 0.0);
            camera.position = camera.position + up;
            cout << "UP" << endl; } 
        else if (event.key.keysym.sym == SDLK_DOWN) {
            glm::vec3 down = glm::vec3(0.0, 0.1, 0.0);
            camera.position = camera.position + down;
            cout << "DOWN" << endl; }
        else if (event.key.keysym.sym == SDLK_p) {
			if (model != "cornell-box") {
				model = "cornell-box";
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0,0.0,1.0); camera.focalLength = 1.0;
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box.obj", 0.17);
			}
            renderMode = "pointCloud";
            cout << "Point Cloud" << endl; }
        else if (event.key.keysym.sym == SDLK_w) {
			if (model != "cornell-box") {
				model = "cornell-box";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0,0.0,1.0); camera.focalLength = 1.0;
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box.obj", 0.17);
			}
            renderMode = "wireFrame";
            cout << "wireFrame" << endl; }
        else if (event.key.keysym.sym == SDLK_r) {
			if (model != "cornell-box") {
				model = "cornell-box";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0,0.0,1.0); camera.focalLength = 1.0;
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box.obj", 0.17);
			}
            renderMode = "rasterise";
            cout << "Rasterise" << endl; }
        else if (event.key.keysym.sym == SDLK_t) {
			if (model != "cornell-box") {
				model = "cornell-box";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0,0.0,1.0); camera.focalLength = 1.0;
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box.obj", 0.17);
			}
            renderMode = "raytrace";
            cout << "Ray Trace" << endl;
        }
		else if (event.key.keysym.sym == SDLK_g) {
			if (model != "sphere") {
				model = "sphere";
				clearAll();
				lightPosition = glm::vec3(0.2, 0.0, 0.5); camera.position = glm::vec3(0.0,0.0,4.0); camera.focalLength = 2.0;
				parseOBJ("build/models/sphere.obj", 0.17);
			}
			renderMode = "gouraud";
			cout << "Gouraud Shading" << endl;
		}
		else if (event.key.keysym.sym == SDLK_h) {
			if (model != "sphere") {
				model = "sphere";
				clearAll();
				lightPosition = glm::vec3(0.2, 0.0, 0.5); camera.position = glm::vec3(0.0,0.0,4.0); camera.focalLength = 2.0;
				parseOBJ("build/models/sphere.obj", 0.17);
			}
			renderMode = "phong";
			cout << "Phong Shading" << endl;
		}
		else if (event.key.keysym.sym == SDLK_m) {
			if (model != "cornell-box-reflection") {
				model = "cornell-box-reflection";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0,0.0,1.0); camera.focalLength = 1.0;
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box-reflection.obj", 0.17);
			}
			renderMode = "reflection";
			cout << "Reflection" << endl;
		}
        draw(window);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
