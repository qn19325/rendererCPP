#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <unordered_map>
#include <RayTriangleIntersection.h>
#include <math.h>

#define WIDTH 320*3
#define HEIGHT 240*3

using namespace std;

int imageNumber = 1550;

string model;
vector<glm::vec3> vertices;
vector<glm::vec3> normals;
vector<ModelTriangle> faces;
vector<CanvasTriangle> canvasTriangles;
unordered_map<string, Colour> colourPalette;
vector<Colour> triangleColours; 
TextureMap texture;
float depthBuffer[WIDTH][HEIGHT];
string renderMode;
glm::vec3 lightPosition;
vector<glm::vec3> lightCluster;
bool orbitMode = false;
struct camera {
    glm::vec3 position;
    float focalLength;
	glm::vec3 right = glm::vec3(1,0,0);
	glm::vec3 up = glm::vec3(0,1,0);
	glm::vec3 forward = glm::vec3(0,0,1);
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
void clearDepthBuffer() {
    for(int i=0; i<WIDTH; i++){
		for(int j=0; j<HEIGHT; j++){
			depthBuffer[i][j] = 0;
		}
	}
}
void lookAt() {
	glm::vec3 forward = glm::normalize(camera.position);
	glm::vec3 vertical(0, 1, 0);
	glm::vec3 right = glm::normalize(glm::cross(vertical, forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));
	glm::mat3 newOrientation(right, up, forward);
	camera.orientation = newOrientation;
}
glm::vec3 getVertexRelativeToCamera(glm::vec3 vertexPosition) {
	glm::vec3 cameraToVertex = vertexPosition - camera.position;
    glm::vec3 adjustedVector = cameraToVertex * camera.orientation;
	return adjustedVector;
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
CanvasPoint getCanvasIntersectionPointTex(glm::vec3 vertexPosition, DrawingWindow &window) {
	glm::vec3 relativePos = getVertexRelativeToCamera(vertexPosition);
	float x = relativePos.x;
	float y = relativePos.y;
	float depth = relativePos.z;
	float u = 400 * (camera.focalLength * (x / depth)) + (window.width / 2);
	float v = 400 * (camera.focalLength * (y / depth)) + (window.height / 2);
	u = WIDTH - u;
	if(depth != 0.0) depth = abs(1 / depth);
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
vector<CanvasPoint> interpolationWithSteps(CanvasPoint to, CanvasPoint from, int steps) {
	vector<CanvasPoint> line;
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
	float xStepSize = xDiff/steps;
	float yStepSize = yDiff/steps;
	float zStepSize = zDiff/steps;
	for (float i=0.0; i<steps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		float z = from.depth + (zStepSize*i);
		line.push_back(CanvasPoint(x, y, z));
	}
	return line;
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
	if (orbitMode == true) {
		glm::mat3 rotation = yRotationMatrix(1.0);
		camera.position = rotation * camera.position; 
		lookAt();
	} else {
		camera.position = glm::vec3(0.0,0.0,3.0); 
		camera.focalLength = 2.0;
		glm::vec3 right = glm::vec3(1,0,0);
		glm::vec3 up = glm::vec3(0,1,0);
		glm::vec3 forward = glm::vec3(0,0,1);
		camera.orientation = glm::mat3(right, up, forward);
	}
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
	// imageNumber += 1;
	// string start;
	// if (imageNumber < 10) {
	// 	start = "images/0000" + to_string(imageNumber) + ".ppm";
	// } else if (imageNumber < 100) {
	// 	start = "images/000" + to_string(imageNumber) + ".ppm";
	// } else if (imageNumber < 1000) {
	// 	start = "images/00" + to_string(imageNumber) + ".ppm";
	// }
	// window.savePPM(start);
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
CanvasTriangle sortTriangleTex(CanvasTriangle triangle) {
	if(triangle.v0().y > triangle.v1().y) {
		swap(triangle.v1(), triangle.v0());
	}
	if(triangle.v0().y > triangle.v2().y) {
		swap(triangle.v2(), triangle.v0());
	}
	if(triangle.v1().y > triangle.v2().y) {
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
CanvasPoint getPoint(CanvasPoint top, CanvasPoint middle, CanvasPoint bottom) {
	float proportion = (bottom.x - top.x) / (bottom.y - top.y);
	float yCoord =  middle.y;
	float xCoord = top.x + (middle.y - top.y) * proportion;
	
	float depth = FLT_MAX;
	float distanceTopBottom = glm::length(glm::vec2(bottom.x - top.x,bottom.y - top.y));
	float distanceTopPoint = glm::length(glm::vec2((xCoord - top.x),(yCoord - top.y)));
	float depthTopBottom = bottom.depth - top.depth;
	float depthTopPoint = (depthTopBottom*distanceTopPoint)/distanceTopBottom ;
	depth = depthTopPoint + top.depth;
	return CanvasPoint(xCoord, yCoord, depth);
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

void filled(CanvasTriangle tri, DrawingWindow &window, Colour col, int side){
	vector<CanvasPoint> side1 = interpolateLine(tri.v0(), tri.v1(), side);
	vector<CanvasPoint> side2 = interpolateLine(tri.v0(), tri.v2(), side);
	for(int i = 0; i< side1.size(); i++){
		lineDrawingDepth(window, side1[i], side2[i], col);
	}
	strokedDepth(window, tri, col);
}
void fillTopTriangleTex(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture) {
	vector<CanvasPoint> side1 = interpolationWithSteps(triangle.v0(), triangle.v1(), (int)(triangle.v1().y- triangle.v0().y)+1);
	vector<CanvasPoint> side2 = interpolationWithSteps(triangle.v0(), triangle.v2(), (int)(triangle.v2().y - triangle.v0().y)+1);
	CanvasPoint v0Tex(triangle.v0().texturePoint.x, triangle.v0().texturePoint.y);
	CanvasPoint v1Tex(triangle.v1().texturePoint.x, triangle.v1().texturePoint.y);
	CanvasPoint v2Tex(triangle.v2().texturePoint.x, triangle.v2().texturePoint.y);
	vector<CanvasPoint> side1Tex = interpolationWithSteps(v0Tex, v1Tex, (int)(triangle.v1().y- triangle.v0().y)+1);
	vector<CanvasPoint> side2Tex = interpolationWithSteps(v0Tex, v2Tex, (int)(triangle.v2().y - triangle.v0().y)+1);

	for(int i = 0; i < side1Tex.size(); i++) {
		float startXTex = side1Tex[i].x; float startYTex = side1Tex[i].y; 
		float endXTex = side2Tex[i].x; float endYTex = side2Tex[i].y; 
		float startX = side1[i].x; float startY = side1[i].y; float startZ = side1[i].depth;
		float endX = side2[i].x; float endY = side2[i].y; float endZ = side2[i].depth;
		
		float xTDiff = endXTex - startXTex; float yTDiff = endYTex - startYTex;
		float xDiff = endX - startX; float yDiff = endY - startY; float zDiff = endZ - startZ;

		float numberOfSteps = max(abs(xDiff), abs(yDiff));
		
		float xStepTex = xTDiff / numberOfSteps; float yStepTex = yTDiff / numberOfSteps;
		float xStep = xDiff / numberOfSteps; float yStep = yDiff / numberOfSteps; float zStep = zDiff / numberOfSteps;

		for(float i = 0; i < numberOfSteps; i++) {
			float xPixelT = startXTex + xStepTex * i; float yPixelT = startYTex + yStepTex * i;
			float xPixel = startX + xStep * i; float yPixel = startY + yStep * i; float zPixel = startZ + zStep * i; 
			uint32_t colour = texture.pixels[(floor(xPixelT) + (texture.width) * floor(yPixelT))];
			int xCoord = floor(xPixel);
			int yCoord = floor(yPixel);
			if(zPixel == 0 || zPixel > (depthBuffer[xCoord][yCoord])) {
				depthBuffer[xCoord][yCoord] = 1/zPixel;
				window.setPixelColour(xCoord, yCoord, colour);
			}
		}
	}
}

void fillBottomTriangleTex(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture) {
	vector<CanvasPoint> side1 = interpolationWithSteps(triangle.v2(), triangle.v1(), (int)(triangle.v2().y - triangle.v1().y)+1);
	vector<CanvasPoint> side2 = interpolationWithSteps(triangle.v2(), triangle.v0(), (int)(triangle.v2().y - triangle.v0().y)+1);
	CanvasPoint v0Tex(triangle.v0().texturePoint.x, triangle.v0().texturePoint.y);
	CanvasPoint v1Tex(triangle.v1().texturePoint.x, triangle.v1().texturePoint.y);
	CanvasPoint v2Tex(triangle.v2().texturePoint.x, triangle.v2().texturePoint.y);
	vector<CanvasPoint> side1Tex = interpolationWithSteps(v2Tex, v1Tex, (int)(triangle.v2().y- triangle.v1().y)+1);
	vector<CanvasPoint> side2Tex = interpolationWithSteps(v2Tex, v0Tex, (int)(triangle.v2().y - triangle.v0().y)+1);

	for(int i = 0; i < side1Tex.size(); i++) {
		float startXTex = side1Tex[i].x; float startYTex = side1Tex[i].y; 
		float endXTex = side2Tex[i].x; float endYTex = side2Tex[i].y; 
		float startX = side1[i].x; float startY = side1[i].y; float startZ = side1[i].depth;
		float endX = side2[i].x; float endY = side2[i].y; float endZ = side2[i].depth;
		
		float xTDiff = endXTex - startXTex; float yTDiff = endYTex - startYTex;
		float xDiff = endX - startX; float yDiff = endY - startY; float zDiff = endZ - startZ;

		float numberOfSteps = max(abs(xDiff), abs(yDiff));
		
		float xStepTex = xTDiff / numberOfSteps; float yStepTex = yTDiff / numberOfSteps;
		float xStep = xDiff / numberOfSteps; float yStep = yDiff / numberOfSteps; float zStep = zDiff / numberOfSteps;

		for(float i = 0; i < numberOfSteps; i++) {
			float xPixelT = startXTex + xStepTex * i; float yPixelT = startYTex + yStepTex * i;
			float xPixel = startX + xStep * i; float yPixel = startY + yStep * i; float zPixel = startZ + zStep * i; 
			uint32_t colour = texture.pixels[(floor(xPixelT) + (texture.width) * floor(yPixelT))];
			int xCoord = floor(xPixel);
			int yCoord = floor(yPixel);
			if(zPixel == 0 || zPixel > (depthBuffer[xCoord][yCoord])) {
				depthBuffer[xCoord][yCoord] = 1/zPixel;
				window.setPixelColour(xCoord, yCoord, colour);
			}	
		}

	}	
}
void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle unsortedTriangle, TextureMap texture) {
	CanvasTriangle triangle = sortTriangleTex(unsortedTriangle);
	CanvasPoint point = getPoint(triangle.v0(), triangle.v1(), triangle.v2());
	CanvasTriangle topTri = CanvasTriangle(triangle.v0(), triangle.v1(), point);
	CanvasTriangle bottomTri = CanvasTriangle(point, triangle.v1(),triangle.v2());
	fillTopTriangleTex(window, topTri, texture);
	fillBottomTriangleTex(window, bottomTri, texture);	
	
}
void rasterise(DrawingWindow &window){
    window.clearPixels();
    clearDepthBuffer();
    canvasTriangles.clear();
    if (canvasTriangles.size() == 0) {
        vector<CanvasPoint> canvasPoints;
        for(int i=0; i<faces.size(); i++) {
            for(int j=0; j<3; j++){
				CanvasPoint temp;
				if (faces[i].colour.name == "Cobbles") temp = getCanvasIntersectionPointTex(faces[i].vertices[j], window);
                else temp = getCanvasIntersectionPoint(faces[i].vertices[j], window);
                canvasPoints.push_back(temp);
            }
			Colour col = faces[i].colour;
			col.name = faces[i].colour.name;
            triangleColours.push_back(col);
        }
        for(int i=0; i<canvasPoints.size(); i=i+3){
            CanvasTriangle temp = CanvasTriangle(canvasPoints[i], canvasPoints[i+1], canvasPoints[i+2]);
            canvasTriangles.push_back(temp);
        }
    }
	if (orbitMode == true) {
		glm::mat3 rotation = yRotationMatrix(1.0);
		camera.position = rotation * camera.position; 
		lookAt();
	} else {
		camera.position = glm::vec3(0.0,0.0,2.0); 
		camera.focalLength = 2.0;
		glm::vec3 right = glm::vec3(1,0,0);
		glm::vec3 up = glm::vec3(0,1,0);
		glm::vec3 forward = glm::vec3(0,0,1);
		camera.orientation = glm::mat3(right, up, forward);
	}
    for (int i=0; i<canvasTriangles.size(); i++) {
		ModelTriangle triangle = faces[i];
        Colour col = triangleColours[i];
        CanvasTriangle canvasTriangle = canvasTriangles[i];
		if (triangleColours[i].name == "Cobbles") {	
			canvasTriangle.v0().texturePoint = TexturePoint(triangle.texturePoints[0].x * texture.width, texture.height - triangle.texturePoints[0].y * texture.height);
			canvasTriangle.v1().texturePoint = TexturePoint(triangle.texturePoints[1].x * texture.width, texture.height - triangle.texturePoints[1].y * texture.height);
			canvasTriangle.v2().texturePoint = TexturePoint(triangle.texturePoints[2].x * texture.width, texture.height - triangle.texturePoints[2].y * texture.height);
			drawTexturedTriangle(window, canvasTriangle, texture);
		} else {
			CanvasTriangle tri = sortTriangle(canvasTriangle);
        	CanvasPoint newp = getMiddlePoint(tri);
			CanvasTriangle topTri = CanvasTriangle(tri.v0(), tri.v1(), newp);
			CanvasTriangle bottomTri = CanvasTriangle(tri.v2(), tri.v1(), newp);
			filled(topTri, window, col, 0);
			filled(bottomTri, window, col, 1);
		}
    }
	// int i = 0;
	// while (i < 75) {
	// 	imageNumber += 1;
	// 	string start;
	// 	if (imageNumber < 10) {
	// 		start = "images/0000" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 100) {
	// 		start = "images/000" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 1000) {
	// 		start = "images/00" + to_string(imageNumber) + ".ppm";
	// 	}
	// 	window.savePPM(start);
	// 	i++;
	// }
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

bool rayCanReachLight(glm::vec3 point, int index, glm::vec3 lightPosition) {
	glm::vec3 rayDirection = lightPosition - point;
	RayTriangleIntersection intersection = getClosestIntersectionRayTrace(rayDirection, point, index);
	if (intersection.triangleIndex != -1) {
		glm::vec3 distanceToIntersection = intersection.intersectionPoint - point;
		glm::vec3 distanceToLight = lightPosition - point ;
		if (glm::length(distanceToIntersection) > glm::length(distanceToLight)) return true;
		else return false;
	}
	return true;
}

float proximity(glm::vec3 lightRay) {
	float proximity = 1 / (2 * M_PI * glm::pow(glm::length(lightRay), 2));
	if (proximity > 1) proximity = 1;
	return proximity;
}
float angleOfIncidence(RayTriangleIntersection intersection, glm::vec3 lightRay) {
	float angleOfIncidence = glm::dot(glm::normalize(intersection.intersectedTriangle.faceNormal), glm::normalize(lightRay));
	if (angleOfIncidence < 0) angleOfIncidence = 0;
	return angleOfIncidence;
}
float specular(RayTriangleIntersection intersection, int n) {
	glm::vec3 lightToPoint = glm::normalize(intersection.intersectionPoint - lightPosition);
	glm::vec3 normalizedNormal = glm::normalize(intersection.intersectedTriangle.faceNormal);
	glm::vec3 rayReflection = glm::normalize(lightToPoint - normalizedNormal * 2.0f * glm::dot(lightToPoint, normalizedNormal)); 
	glm::vec3 toCamera = glm::normalize(camera.position - intersection.intersectionPoint);
	float specularLight = glm::dot(toCamera,rayReflection);
	if (specularLight < 0) specularLight = 0;
	specularLight = pow(specularLight, n);
	return specularLight;
}

void raytrace(DrawingWindow& window) {
	window.clearPixels();
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			float xRayDirection = x - WIDTH / 2.0;
			float yRayDirection = y - HEIGHT / 2.0;
			xRayDirection = xRayDirection / 400;
			yRayDirection = -yRayDirection / 400;
			glm::vec3 imagePlanePoint(xRayDirection, yRayDirection, camera.position.z - camera.focalLength);
			glm::vec3 rayDirection = imagePlanePoint - camera.position;
			rayDirection = glm::normalize(rayDirection);
			RayTriangleIntersection intersection = getClosestIntersectionRayTrace(rayDirection, camera.position, -1);
			if (intersection.intersectedTriangle.reflective == true) {
				glm::vec3 reflectionVector = rayDirection - (2.0f * intersection.intersectedTriangle.faceNormal * glm::dot(rayDirection, intersection.intersectedTriangle.faceNormal));
				RayTriangleIntersection reflectionIntersection = getClosestIntersectionRayTrace(glm::normalize(reflectionVector), intersection.intersectionPoint + 0.0001f * intersection.intersectedTriangle.faceNormal, -1);
				intersection = reflectionIntersection;
			}
			if (intersection.triangleIndex != -1) {
				Colour colour;
                glm::vec3 lightRay = glm::vec3(lightPosition - intersection.intersectionPoint);
				
                float prox = proximity(lightRay);
				float aoi = angleOfIncidence(intersection, lightRay);				
				float spec = specular(intersection, 10);
				
				float brightness = (prox + aoi + spec) / 3;
				if ((brightness+0.2) < 1) brightness += 0.2;
				else brightness = 1.0;
				if (rayCanReachLight(intersection.intersectionPoint,intersection.triangleIndex, lightPosition)) {
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
	while (imageNumber < 1350) {
		imageNumber += 1;
		string start;
		if (imageNumber < 10) {
			start = "images/0000" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 100) {
			start = "images/000" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 1000) {
			start = "images/00" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 10000) {
			start = "images/0" + to_string(imageNumber) + ".ppm";
		}
		window.savePPM(start);
	}
}
void addLights() {
	lightCluster.clear();
	for (float i=0.0; i<0.08; i=i+0.01) {
		// create a cluster of light points around the original light
		lightCluster.push_back(glm::vec3(lightPosition.x + i, lightPosition.y, lightPosition.z));
		lightCluster.push_back(glm::vec3(lightPosition.x - i, lightPosition.y, lightPosition.z));
		lightCluster.push_back(glm::vec3(lightPosition.x, lightPosition.y, lightPosition.z + i));
		lightCluster.push_back(glm::vec3(lightPosition.x, lightPosition.y, lightPosition.z - i));
		lightCluster.push_back(glm::vec3(lightPosition.x + i, lightPosition.y, lightPosition.z + i));
		lightCluster.push_back(glm::vec3(lightPosition.x - i, lightPosition.y, lightPosition.z - i));
		lightCluster.push_back(glm::vec3(lightPosition.x - i, lightPosition.y, lightPosition.z + i));
		lightCluster.push_back(glm::vec3(lightPosition.x + i, lightPosition.y, lightPosition.z - i));
	}
}
void softShadow(DrawingWindow& window) {
	window.clearPixels();
	addLights();
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			float xRayDirection = x - WIDTH / 2.0;
			float yRayDirection = y - HEIGHT / 2.0;
			xRayDirection = xRayDirection / 400;
			yRayDirection = -yRayDirection / 400;
			glm::vec3 imagePlanePoint(xRayDirection, yRayDirection, camera.position.z - camera.focalLength);
			glm::vec3 rayDirection = imagePlanePoint - camera.position;
			rayDirection = glm::normalize(rayDirection);
			RayTriangleIntersection intersection = getClosestIntersectionRayTrace(rayDirection, camera.position, -1);
			if (intersection.triangleIndex != -1) {
				Colour colour;
                glm::vec3 lightRay = glm::vec3(lightPosition - intersection.intersectionPoint);
				
				float prox = proximity(lightRay);
				float aoi = angleOfIncidence(intersection, lightRay);				
				float spec = specular(intersection, 10);
				
				float brightness = (prox + aoi + spec) / 3;
				float numberHit = 0;
				for (int i=0; i<lightCluster.size(); i++) {
					if (rayCanReachLight(intersection.intersectionPoint, intersection.triangleIndex, lightCluster[i])) {
						numberHit++;
					}
				}
				float intensity = numberHit / lightCluster.size();
				if (intensity < 0.2) intensity =0.2;
				colour = intersection.intersectedTriangle.colour;
				uint32_t col = (255 << 24) + (int(colour.red*intensity) << 16) + (int(colour.green*intensity) << 8) + int(colour.blue*intensity);
				window.setPixelColour(x, y, col);
			}
		}
	}
	// while (imageNumber < 1450) {
	// 	imageNumber += 1;
	// 	string start;
	// 	if (imageNumber < 10) {
	// 		start = "images/0000" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 100) {
	// 		start = "images/000" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 1000) {
	// 		start = "images/00" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 10000) {
	// 		start = "images/0" + to_string(imageNumber) + ".ppm";
	// 	}
	// 	window.savePPM(start);
	// }
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
					int n = 256;
					glm::vec3 light2Point = glm::normalize(intersection.intersectionPoint - lightPosition);
					normal = glm::normalize(normal);
					glm::vec3 reflection = glm::normalize((light2Point - 2.0f*normal*(glm::dot(light2Point, normal))));
					glm::vec3 toCamera = glm::normalize((camera.position - intersection.intersectionPoint));
					float spec = glm::dot(toCamera, reflection);
					if (spec < 0) spec = 0;
					spec = pow(spec, n);
					float brightness = (aoi + spec)/2;
					face.vertexBrightness[i] = brightness;
				}
				float brightness = (face.vertexBrightness[0]*intersection.baryCoords[0]) + (face.vertexBrightness[1]*intersection.baryCoords[1]) + (face.vertexBrightness[2]*intersection.baryCoords[2]);
                Colour colour = Colour(255,0,0);
				uint32_t col = (255 << 24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + int(colour.blue*brightness);
				window.setPixelColour(x, y, col);
			}
		}
	}
	// while (imageNumber < 1550) {
	// 	imageNumber += 1;
	// 	string start;
	// 	if (imageNumber < 10) {
	// 		start = "images/0000" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 100) {
	// 		start = "images/000" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 1000) {
	// 		start = "images/00" + to_string(imageNumber) + ".ppm";
	// 	} else if (imageNumber < 10000) {
	// 		start = "images/0" + to_string(imageNumber) + ".ppm";
	// 	}
	// 	window.savePPM(start);
	// }
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
				glm::vec3 pixelNormal1 = glm::normalize(intersection.intersectedTriangle.vertexNormals[0])*intersection.baryCoords[0];
				glm::vec3 pixelNormal2 = glm::normalize(intersection.intersectedTriangle.vertexNormals[1])*intersection.baryCoords[1];
				glm::vec3 pixelNormal3 = glm::normalize(intersection.intersectedTriangle.vertexNormals[2])*intersection.baryCoords[2];
				glm::vec3 pixelNormal = glm::normalize(pixelNormal1 + pixelNormal2 + pixelNormal3);

				glm::vec3 vectorFromLight = glm::normalize(lightPosition - intersection.intersectionPoint);
				float aoi = glm::dot(pixelNormal, vectorFromLight);
				if(aoi<0) aoi = 0;

				glm::vec3 light2Point = glm::normalize(intersection.intersectionPoint - lightPosition);
				glm::vec3 reflection = glm::normalize((light2Point - 2.0f*pixelNormal*(glm::dot(light2Point, pixelNormal))));
				glm::vec3 pointToCamera = glm::normalize((camera.position - intersection.intersectionPoint));
				float specular = glm::dot(pointToCamera, reflection);
				if (specular < 0) specular = 0;
				specular = pow(specular, 256);
				
				float brightness = (specular + aoi);
				if (brightness<0) brightness = 0;
				if (brightness>1) brightness = 1;
				Colour colour = Colour(255,0,0);
				uint32_t col = (255<<24) + (int(colour.red*brightness) << 16) + (int(colour.green*brightness) << 8) + (int(colour.blue*brightness));
				window.setPixelColour(x, y, col);
			}

		}
	}
	while (imageNumber < 1650) {
		imageNumber += 1;
		string start;
		if (imageNumber < 10) {
			start = "images/0000" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 100) {
			start = "images/000" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 1000) {
			start = "images/00" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 10000) {
			start = "images/0" + to_string(imageNumber) + ".ppm";
		}
		window.savePPM(start);
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
	vector<TexturePoint> texturePoints;
	bool reflective = false;
	fileIn.open(filePath);
	if (fileIn.is_open()) {
		while (getline(fileIn, line)) {
			if (line[0] == 'v') {
				vector<string> lineParts = split(line, ' ');
				if (line[1] == 'n') {
					glm::vec3 normal = glm::vec3(stof(lineParts[1]),stof(lineParts[2]),stof(lineParts[3]));
					normals.push_back(normal);
				} else if (line[1] == 't') {
					TexturePoint tp = TexturePoint(stof(lineParts[1]),stof(lineParts[2]));
					texturePoints.push_back(tp);
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
				colour.name = currentColour;
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
				if (colour.name == "Cobbles") {
					face.texturePoints[0] = texturePoints[stoi(splitLineParts1[1])-1];
					face.texturePoints[1] = texturePoints[stoi(splitLineParts2[1])-1];
					face.texturePoints[2] = texturePoints[stoi(splitLineParts3[1])-1];
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
			} else if(line[0] == 'm') {
				vector<string> lineParts = split(line, ' ');
				texture = TextureMap("build/models/" + lineParts[1]);
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
	else if (renderMode == "rasteriseTexture") rasterise(window);
	else if (renderMode == "gouraud") gouraud(window);
	else if (renderMode == "phong") phong(window);
	else if (renderMode == "reflection") raytrace(window);
	else if (renderMode == "softShadows") softShadow(window);
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
			lightPosition = lightPosition + down;
            cout << "DOWN" << endl; }
		else if (event.key.keysym.sym == SDLK_x) {
            glm::mat3 rotationMatrix = xRotationMatrix(1.0);
            camera.position = camera.position * rotationMatrix;
            cout << "ROTATE X" << endl;
        } else if (event.key.keysym.sym == SDLK_y) {
            glm::mat3 rotationMatrix = yRotationMatrix(1.0);
            camera.position = camera.position * rotationMatrix;
            cout << "ROTATE Y" << endl;
        } else if (event.key.keysym.sym == SDLK_EQUALS) {
            glm::mat3 rotationMatrix = yRotationMatrix(1.0);
            camera.orientation = rotationMatrix * camera.orientation;
        } else if (event.key.keysym.sym == SDLK_MINUS) {
            glm::mat3 rotationMatrix = xRotationMatrix(1.0);
            camera.orientation = rotationMatrix * camera.orientation;
        } else if (event.key.keysym.sym == SDLK_o) {
			if (orbitMode == true) orbitMode = false;
			else orbitMode = true;
		}
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
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box.obj", 0.17);
			}
            renderMode = "rasterise";
            cout << "Rasterise" << endl; }
		else if (event.key.keysym.sym == SDLK_c) {
			if (model != "textured-cornell-box") {
				model = "textured-cornell-box";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0,0.0,1.0); camera.focalLength = 1.0;
				parseMTL("build/models/textured-cornell-box.mtl");
				parseOBJ("build/models/textured-cornell-box.obj", 0.17);
			}
            renderMode = "rasteriseTexture";
            cout << "Rasterise" << endl; }
        else if (event.key.keysym.sym == SDLK_t) {
			if (model != "cornell-box") {
				model = "cornell-box";
				clearAll();
				window.clearPixels();
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
				lightPosition = glm::vec3(0.0, 0.9, 0.9); camera.position = glm::vec3(0.0,0.0,4.0); camera.focalLength = 2.0;
				parseOBJ("build/models/sphere.obj", 0.07);
			}
			renderMode = "gouraud";
			cout << "Gouraud Shading" << endl;
		}
		else if (event.key.keysym.sym == SDLK_h) {
			if (model != "sphere") {
				model = "sphere";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.9, 0.9); camera.position = glm::vec3(0.0,0.0,4.0); camera.focalLength = 2.0;
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
		else if (event.key.keysym.sym == SDLK_s) {
			if (model != "cornell-box") {
				model = "cornell-box";
				clearAll();
				lightPosition = glm::vec3(0.0, 0.4, 0.0); camera.position = glm::vec3(0.0, 0.0, 1.0); camera.focalLength = 1.0;
				parseMTL("build/models/cornell-box.mtl");
				parseOBJ("build/models/cornell-box.obj", 0.17);
			}
			renderMode = "softShadows";
			cout << "Soft Shadows" << endl;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		imageNumber += 1;
		string start;
		if (imageNumber < 10) {
			start = "images/0000" + to_string(imageNumber) + ".ppm";
		} else if (imageNumber < 100) {
			start = "images/000" + to_string(imageNumber) + ".ppm";
		}
		window.savePPM(start);
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
