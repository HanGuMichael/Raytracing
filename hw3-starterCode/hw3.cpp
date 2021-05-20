/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Han Gu
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>


#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100
#define PI 3.1415926535

#define RECURSIVE 0

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];
int render[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Point {
    double x;
    double y;
    double z;
    Point() {
        x = 0;
        y = 0;
        z = 0;
    }
    Point(double a,double b,double c) {
        x = a;
        y = b;
        z = c;
    }
    

    double length() {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        double l = length();
        x = x / l;
        y = y / l;
        z = z / l;
    }

    double dot(Point a) {
        return x * a.x + y * a.y + z * a.z;
    }

    Point cross(Point b) {
        Point ret;
        ret.x = y * b.z - z * b.y;
        ret.y = z * b.x - x * b.z;
        ret.z = x * b.y - y * b.x;
        return ret;
    }

    Point operator*(double a) {
        return Point(x * a, y * a, z * a);
    }

    Point operator+(Point p1) {
        return Point(x + p1.x, y + p1.y, z + p1.z);
    }

    Point operator-(Point p1) {
        return Point(x - p1.x, y - p1.y, z - p1.z);
    }
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//intersection funcitons
double sphereIntersection(Point p0, Point pd, int sphere) {
    double xc = spheres[sphere].position[0];
    double yc = spheres[sphere].position[1];
    double zc = spheres[sphere].position[2];
    double r = spheres[sphere].radius;

    double b = 2 * (pd.x * (p0.x - xc) + pd.y * (p0.y - yc) + pd.z * (p0.z - zc));
    double c = (p0.x - xc) * (p0.x - xc) + (p0.y - yc) * (p0.y - yc) + (p0.z - zc) * (p0.z - zc) - r * r;

    if (b * b > 4 * c) {
        double t0 = (-b - (sqrt(b * b - 4 * c))) / 2;
        double t1 = (-b + (sqrt(b * b - 4 * c))) / 2;
        if (t0 > 0.0001) {
            return t0;
        }
        else if (t1 > 0.0001) {
            return t1;
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }

}

double area(Point a, Point b, Point c) {
    return (double)1.0 / 2 * ((b - a).cross(c - a)).length();
}

double triangleIntersection(Point p0, Point pd, int triangle, Point & np, Point & diffuse, Point & specular, double &shiness) { //np is normal at intersection point
    Point a = Point(triangles[triangle].v[0].position[0], triangles[triangle].v[0].position[1], triangles[triangle].v[0].position[2]);
    Point b = Point(triangles[triangle].v[1].position[0], triangles[triangle].v[1].position[1], triangles[triangle].v[1].position[2]);
    Point c = Point(triangles[triangle].v[2].position[0], triangles[triangle].v[2].position[1], triangles[triangle].v[2].position[2]);

    Point normal = (a - b).cross(b - c);
    normal.normalize();

    double d = -(normal.x * a.x) - (normal.y * a.y) - (normal.z * a.z);
    double t = -(normal.x * p0.x + normal.y * p0.y + normal.z * p0.z + d) / (normal.x * pd.x + normal.y * pd.y + normal.z * pd.z);
    if (t > 0.0001) {
        //intersect the plane,check if inside triangle
        Point C = Point(p0.x + pd.x * t, p0.y + pd.y * t, p0.z + pd.z * t);
        double areaTotal = area(a, b, c);
        double alpha = area(C, b, c) / areaTotal;
        double beta = area(a, C, c) / areaTotal;
        double gamma = area(a, b, C) / areaTotal;
        //double gamma = 1 - alpha - beta;
        if (alpha < 0 || alpha>1) {
            return 0;
        }
        else if (beta < 0 || beta>1) {
            return 0;
        }
        else if (gamma < 0 || gamma>1) {
            return 0;
        }
        else if (alpha + beta + gamma > 1.000001 || alpha + beta + gamma < 0.999999) {
            return 0;
        }
        Point n0 = Point(triangles[triangle].v[0].normal[0], triangles[triangle].v[0].normal[1], triangles[triangle].v[0].normal[2]);
        Point n1= Point(triangles[triangle].v[1].normal[0], triangles[triangle].v[1].normal[1], triangles[triangle].v[1].normal[2]);
        Point n2 = Point(triangles[triangle].v[2].normal[0], triangles[triangle].v[2].normal[1], triangles[triangle].v[2].normal[2]);
        Point temp = n0 * alpha + n1 * beta + n2 * gamma;
        temp.normalize();
        np.x = temp.x;
        np.y = temp.y;
        np.z = temp.z;
        diffuse.x = triangles[triangle].v[0].color_diffuse[0] * alpha + triangles[triangle].v[1].color_diffuse[0] * beta + triangles[triangle].v[2].color_diffuse[0] * gamma;
        diffuse.y = triangles[triangle].v[0].color_diffuse[1] * alpha + triangles[triangle].v[1].color_diffuse[1] * beta + triangles[triangle].v[2].color_diffuse[1] * gamma;
        diffuse.z = triangles[triangle].v[0].color_diffuse[2] * alpha + triangles[triangle].v[1].color_diffuse[2] * beta + triangles[triangle].v[2].color_diffuse[2] * gamma;
        specular.x = triangles[triangle].v[0].color_specular[0] * alpha + triangles[triangle].v[1].color_specular[0] * beta + triangles[triangle].v[2].color_specular[0] * gamma;
        specular.y = triangles[triangle].v[0].color_specular[1] * alpha + triangles[triangle].v[1].color_specular[1] * beta + triangles[triangle].v[2].color_specular[1] * gamma;
        specular.z = triangles[triangle].v[0].color_specular[2] * alpha + triangles[triangle].v[1].color_specular[2] * beta + triangles[triangle].v[2].color_specular[2] * gamma;
        shiness = triangles[triangle].v[0].shininess * alpha + triangles[triangle].v[1].shininess * beta + triangles[triangle].v[2].shininess * gamma;
        return t;
    }
    else {
        return 0;
    }

}

bool isInShadow(Point p0,int light) {
    Point pd = Point(lights[light].position[0], lights[light].position[1], lights[light].position[2]) - p0;
    double dist = pd.length();
    pd.normalize();
    for (int i = 0; i < num_spheres; i++) {
        double currT = sphereIntersection(p0, pd, i);
        if (currT > 0.001 && currT < dist) { return true; }
    }
    //TODO: loop for triangles
    for (int i = 0; i < num_triangles; i++) {
        Point p1;
        Point p2;
        Point p3;
        double p5 = 0;
        double currT = triangleIntersection(p0, pd, i,p1,p2,p3,p5);
        if (currT > 0.001 && currT < dist) { return true; }
    }

    return false;
}

//double obstacleDistFromIntersection(Point p0, int light , int & obstacleNum, Point & shadowRay) {
//    Point pd = Point(lights[light].position[0], lights[light].position[1], lights[light].position[2]) - p0;
//    pd.normalize();
//    for (int i = 0; i < num_spheres; i++) {
//        double currT = sphereIntersection(p0, pd, i);
//        if (currT > 0.001) { return currT; }
//    }
//    //TODO: loop for triangles
//    for (int i = 0; i < num_triangles; i++) {
//        Point p1;
//        Point p2;
//        Point p3;
//        double p5 = 0;
//        double currT = triangleIntersection(p0, pd, i, p1, p2, p3, p5);
//        if (currT > 0.001) { 
//            obstacleNum = i;
//            shadowRay = pd;
//            return currT; }
//    }
//}
double singleChannelColor(int light, Point p0, Point normal, int channel, double kd, double ks, double sh);
double singleChannelColor2(int light, Point p0, Point normal, int channel, double kd, double ks, double sh);
double reflectSingleChannelColor(int light, int reflectTime, Point p0, Point pd, int channel) {
    if (reflectTime == 0) {
        return 0;
    }
    double smallestT = 99999999;
    bool found = false;
    int sphereIndex = -1;
    int triangleIndex = -1;
    bool onSphere = false;
    Point triNormal0;
    Point diffuse0;
    Point specular0;
    double shiness0 = 0;
    Point  triangleIntersectNormal;
    Point diffuse;
    Point specular;
    double shiness = 0;

    for (unsigned int i = 0; i < num_spheres; i++) {
        double currT = sphereIntersection(p0, pd, i);
        if (currT > 0.001 && currT < smallestT) {
            smallestT = currT;
            found = true;
            onSphere = true;
            sphereIndex = i;
        }
    }
    for (unsigned int i = 0; i < num_triangles; i++) {
        double currT = triangleIntersection(p0, pd, i, triNormal0, diffuse0, specular0, shiness0);
        if (currT > 0.001 && currT < smallestT) {
            triangleIntersectNormal = triNormal0;
            diffuse = diffuse0;
            specular = specular0;
            shiness = shiness0;

            smallestT = currT;
            found = true;
            onSphere = false;
            triangleIndex = i;
        }
    }
    Point intersect = pd * smallestT;
    double ret = 0;
    double coe = 0.5;
    if (found) {
        //return 0.1;
        if (!isInShadow(intersect, light)) {
            //return color
            if (onSphere) {
                Point normal = intersect - Point(spheres[sphereIndex].position[0], spheres[sphereIndex].position[1], spheres[sphereIndex].position[2]);
                normal.normalize();
                Point d = normal * ((pd*(-1)).dot(normal)) * 2 - (pd * (-1));
                d.normalize();
                return coe * reflectSingleChannelColor(light,reflectTime-1,intersect,d,channel)+ (1 - coe)* singleChannelColor2(light, intersect, normal, channel, spheres[sphereIndex].color_diffuse[channel], spheres[sphereIndex].color_specular[channel], spheres[sphereIndex].shininess);
            }
            else {
                Point d = triangleIntersectNormal * ((pd * (-1)).dot(triangleIntersectNormal)) * 2 - (pd * (-1));
                d.normalize();
                if (channel == 0) {
                    return coe*reflectSingleChannelColor(light, reflectTime - 1, intersect, d, channel) + (1-coe)*singleChannelColor2(light, intersect, triangleIntersectNormal, channel, diffuse.x, specular.x, shiness);
                }
                else if (channel == 1) {
                    return coe*reflectSingleChannelColor(light, reflectTime - 1, intersect, d, channel) + (1 - coe) * singleChannelColor2(light, intersect, triangleIntersectNormal, channel, diffuse.y, specular.y, shiness);
                }
                else {
                    return coe*reflectSingleChannelColor(light, reflectTime - 1, intersect, d, channel) + (1 - coe) * singleChannelColor2(light, intersect, triangleIntersectNormal, channel, diffuse.z, specular.z, shiness);
                }
                
            }
        }
        else {
            if (onSphere) {
                //calculate normal & illumination for sphere surface
                Point normal = intersect - Point(spheres[sphereIndex].position[0], spheres[sphereIndex].position[1], spheres[sphereIndex].position[2]);
                normal.normalize();
                Point d = normal * ((pd * (-1)).dot(normal)) * 2 - (pd * (-1));
                d.normalize();
                return coe * reflectSingleChannelColor(light, reflectTime - 1, intersect, d, channel);
            }
            else {
                //intersect = intersect + triangleIntersectNormal * 0.001;
                Point d = triangleIntersectNormal * ((pd * (-1)).dot(triangleIntersectNormal)) * 2 - (pd * (-1));
                d.normalize();
                return coe * reflectSingleChannelColor(light, reflectTime - 1, intersect, d, channel);
            }
        }
    }
    else {
        return 0;
    }
}

double singleChannelColor(int light, Point p0, Point normal,int channel, double kd, double ks, double sh) {
    double ret;
    Point eyeDir;
    if (light != -1) {
        Point lightDir = Point(lights[light].position[0], lights[light].position[1], lights[light].position[2]) - p0;
        lightDir.normalize();
        Point reflectDir = normal * (lightDir.dot(normal)) * 2 - lightDir;
        reflectDir.normalize();
        eyeDir = Point(0, 0, 0) - p0;
        eyeDir.normalize();
        double d = max(lightDir.dot(normal), 0.0f);
        double s = max(reflectDir.dot(eyeDir), 0.0f);
        //return 100;
        ret = lights[light].color[channel] * (kd * d + ks * pow(s, sh));
        if (ret < 0) {
            printf("negative color\n");
        }
    }
    else {
        ret = 0;
    }
    
    //recursive reflection
#if RECURSIVE
    Point pd = normal * (eyeDir.dot(normal)) * 2 - eyeDir;
    double coe = 0.5;
    ret = (1 - coe) * ret + coe * reflectSingleChannelColor(light, 6, p0, pd, channel);

#endif // DEBUG

    return ret;
}

double singleChannelColor2(int light, Point p0, Point normal, int channel, double kd, double ks, double sh) {
    double ret;
    Point eyeDir;
    if (light != -1) {
        Point lightDir = Point(lights[light].position[0], lights[light].position[1], lights[light].position[2]) - p0;
        lightDir.normalize();
        Point reflectDir = normal * (lightDir.dot(normal)) * 2 - lightDir;
        reflectDir.normalize();
        eyeDir = Point(0, 0, 0) - p0;
        eyeDir.normalize();
        double d = max(lightDir.dot(normal), 0.0f);
        double s = max(reflectDir.dot(eyeDir), 0.0f);
        //return 100;
        ret = lights[light].color[channel] * (kd * d + ks * pow(s, sh));
        if (ret < 0) {
            printf("negative color\n");
        }
    }
    else {
        ret = 0;
    }

    return ret;
}


Point PhongShading(Point p0,Point normal, Point kd, Point ks, double sh, Point eye) {
    Point ret = Point(0,0,0);
    Point eyeDir;
    for (int i = 0; i < num_lights; i++) {
        if (!isInShadow(p0, i)) {
            Point lightDir = Point(lights[i].position[0], lights[i].position[1], lights[i].position[2]) - p0;
            lightDir.normalize();
            Point reflectDir = normal * (lightDir.dot(normal)) * 2 - lightDir;
            reflectDir.normalize();
            eyeDir = eye - p0;
            eyeDir.normalize();
            double d = max(lightDir.dot(normal), 0.0f);
            double s = max(reflectDir.dot(eyeDir), 0.0f);
            ret.x += max(0, lights[i].color[0] * (kd.x * d + ks.x * pow(s, sh)));
            ret.y += max(0, lights[i].color[1] * (kd.y * d + ks.y * pow(s, sh)));
            ret.z += max(0, lights[i].color[2] * (kd.z * d + ks.z * pow(s, sh)));
        }
    }
    ret.x = min(1, ret.x);
    ret.y = min(1, ret.y);
    ret.z = min(1, ret.z);
    return ret;
}

Point raytrace(Point p0, Point pd, int recurse, int x, int y) { // returns the rgb(0,1) value observed from p0 to pd. recurse times
    if (recurse == 0) {
        //printf("raytrace recurse: 0\n");
        //printf("coord: %d %d\n",x,y);

        //return Point(100, 0, 0);
        return Point(0, 0, 0);
    }

    double smallestT = 99999999;
    bool found = false;
    int sphereIndex = -1;
    int triangleIndex = -1;
    bool onSphere = false;
    Point triNormal0;
    Point diffuse0;
    Point specular0;
    double shiness0 = 0;
    Point  triangleIntersectNormal;
    Point diffuse;
    Point specular;
    double shiness = 0;

    for (unsigned int i = 0; i < num_spheres; i++) {
        double currT = sphereIntersection(p0, pd, i);
        if (currT > 0.001 && currT < smallestT) {
            smallestT = currT;
            found = true;
            onSphere = true;
            sphereIndex = i;
        }
    }
    for (unsigned int i = 0; i < num_triangles; i++) {
        double currT = triangleIntersection(p0, pd, i, triNormal0, diffuse0, specular0, shiness0);
        if (currT > 0.001 && currT < smallestT) {
            triangleIntersectNormal = triNormal0;
            diffuse = diffuse0;
            specular = specular0;
            shiness = shiness0;

            smallestT = currT;
            found = true;
            onSphere = false;
            triangleIndex = i;
        }
    }

    if (found) {
        //printf("intersect\n");
        //return Point(1, 0, 0);
        Point normal;
        Point kd;
        Point ks;
        double sh;
        Point intersect = pd * smallestT ;
        if (onSphere) {
            //printf("OnSphere\n");
            normal =  intersect - Point(spheres[sphereIndex].position[0], spheres[sphereIndex].position[1], spheres[sphereIndex].position[2]);
            normal.normalize();
            intersect = normal * 0.001 + intersect;
            kd = Point(spheres[sphereIndex].color_diffuse[0], spheres[sphereIndex].color_diffuse[1], spheres[sphereIndex].color_diffuse[2]);
            ks = Point(spheres[sphereIndex].color_specular[0], spheres[sphereIndex].color_specular[1], spheres[sphereIndex].color_specular[2]);
            sh = spheres[sphereIndex].shininess;
        }
        else {
            normal = triangleIntersectNormal;
            normal.normalize();
            kd = diffuse;
            ks = specular;
            sh = shiness;
        }
        Point newRay = normal * (((pd*(-1)).dot(normal)) * 2) - (pd*(-1));
        newRay.normalize();
        double coe = 0.3;
        
        //printf("recurse>1\n");
        Point reflect = raytrace(intersect, newRay, recurse - 1,x,y);
        Point Phong = PhongShading(intersect, normal, kd, ks, sh, p0);
        if (recurse == 1) {
            return Point(Phong.x , Phong.y ,Phong.z);
        }
        //ks = Point(0.1, 0.1, 0.1);
        //ks = kd;
        //return PhongShading(intersect, normal, kd, ks, sh)  + Point(temp.x * ks.x, temp.y * ks.y, temp.z * ks.z);
        return Point(Phong.x * (1-ks.x), Phong.y * (1 - ks.y), Phong.z * (1 - ks.z)) + Point(reflect.x * ks.x, reflect.y * ks.y, reflect.z * ks.z);
    }
    return Point(0, 0, 0);


}


//MODIFY THIS FUNCTION
void draw_scene()
{   
    int renderProgress0 = 0;
    int renderProgress1 = 0;
    double a = (double)WIDTH / (double)HEIGHT;
    for (unsigned int x = 0; x < WIDTH; x++) {
        for (unsigned int y = 0; y < HEIGHT; y++){
            //render[y][x][0] = 255;
            //render[y][x][1] = 0;
            //render[y][x][2] = 0;
            render[y][x][0] = 255;
            render[y][x][1] = 255;
            render[y][x][2] = 255;
            Point p0 = Point(0, 0, 0);
            double xd;
            double yd;
            xd = (x - ((double)WIDTH / 2))/((double)WIDTH / 2) * a * tan(fov*PI / 180 / 2);
            yd = (y - ((double)HEIGHT / 2))/ ((double)HEIGHT / 2) * tan(fov*PI / 180 / 2);
            
            Point pd = Point(xd, yd, -1) - p0;
            pd.normalize();

#if RECURSIVE
            Point rgb = raytrace(p0, pd, 5, x, y);
#else
            Point rgb = raytrace(p0, pd, 1, x, y);
#endif // 0

            render[y][x][0] = rgb.x * 255;
            render[y][x][1] = rgb.y * 255;
            render[y][x][2] = rgb.z * 255;

            /**/
            //double smallestT = 99999999;
            //bool found = false;
            //int sphereIndex = -1;
            //int triangleIndex = -1;
            //bool onSphere = false;
            //Point triNormal0;
            //Point diffuse0;
            //Point specular0;
            //double shiness0 = 0;
            //Point  triangleIntersectNormal;
            //Point diffuse;
            //Point specular;
            //double shiness = 0;
            //
            //for (unsigned int i = 0; i < num_spheres; i++) {
            //    double currT = sphereIntersection(p0, pd, i);
            //    if (currT > 0.001 && currT < smallestT) {
            //        smallestT = currT;
            //        found = true;
            //        onSphere = true;
            //        sphereIndex = i;
            //    }
            //}
            //for (unsigned int i = 0; i < num_triangles; i++) {
            //    double currT = triangleIntersection(p0, pd, i, triNormal0,diffuse0,specular0,shiness0);
            //    if (currT > 0.001 && currT < smallestT) {
            //        triangleIntersectNormal = triNormal0;
            //        diffuse = diffuse0;
            //        specular = specular0;
            //        shiness = shiness0;

            //        smallestT = currT;
            //        found = true;
            //        onSphere = false;
            //        triangleIndex = i;
            //    }
            //}
            ////TODO: intersection with traingles


            ////closest intersection found 

            //if (found) {
            //    render[y][x][0] = 0;
            //    render[y][x][1] = 0;
            //    render[y][x][2] = 0;
            //    //printf("found!\n");
            //    Point intersect = pd * smallestT ;
            //    //intersect = intersect
            //    for (int i = 0; i < num_lights; i++) {
            //        if (!isInShadow(intersect, i)) {
            //            //printf("not in shadow\n");
            //            if (onSphere) {
            //                //calculate normal & illumination for sphere surface
            //                Point normal = intersect - Point(spheres[sphereIndex].position[0], spheres[sphereIndex].position[1], spheres[sphereIndex].position[2]);
            //                normal.normalize();

            //                render[y][x][0] += max(0,255 * singleChannelColor(i, intersect, normal, 0, spheres[sphereIndex].color_diffuse[0], spheres[sphereIndex].color_specular[0], spheres[sphereIndex].shininess));
            //                render[y][x][1] += max(0, 255 * singleChannelColor(i, intersect, normal, 1, spheres[sphereIndex].color_diffuse[1], spheres[sphereIndex].color_specular[1], spheres[sphereIndex].shininess));
            //                render[y][x][2] += max(0, 255 * singleChannelColor(i, intersect, normal, 2, spheres[sphereIndex].color_diffuse[2], spheres[sphereIndex].color_specular[2], spheres[sphereIndex].shininess));
            //            }
            //            else {
            //                intersect = intersect + triangleIntersectNormal * 0.001;
            //                render[y][x][0] += max(0, 255 * singleChannelColor(i, intersect, triangleIntersectNormal, 0, diffuse.x, specular.x, shiness));
            //                render[y][x][1] += max(0, 255 * singleChannelColor(i, intersect,triangleIntersectNormal, 1, diffuse.y, specular.y, shiness));
            //                render[y][x][2] += max(0, 255 * singleChannelColor(i, intersect, triangleIntersectNormal, 2, diffuse.z, specular.z, shiness));
            //            }
            //        }

                    //else {
                    //    //render[y][x][0] = 255;
                    //    if (!onSphere&&triangleIndex==0) {
                    //        //double currT = triangleIntersection(p0, pd, i, triangleIntersectNormal, diffuse, specular, shiness);
                    //        render[y][x][0] = 255;
                    //        Point shadowRay;
                    //        int obstacle = -1;
                    //        printf("Obstacle dist from intersection: %f  ", obstacleDistFromIntersection(intersect,i,obstacle,shadowRay));
                    //        printf("obstacle index: %d shadowRay: %f %f %f\n", obstacle, shadowRay.x, shadowRay.y, shadowRay.z);

                    //    }
                    //}
                    //if in shadow add nothing to render[y][x]
                //}
                /*if (!onSphere) {
                    render[y][x][0] = 100;
                }*/
                render[y][x][0] += 255 * ambient_light[0];
                render[y][x][1] += 255 * ambient_light[0];
                render[y][x][2] += 255 * ambient_light[0];
                if (render[y][x][0] > 255) {
                    render[y][x][0] = 255;
                }
                if (render[y][x][1] > 255) {
                    render[y][x][1] = 255;
                }
                if (render[y][x][2] > 255) {
                    render[y][x][2] = 255;
                }
            
            renderProgress1 = x * 100 / WIDTH;
            if (renderProgress1 > renderProgress0) {
                printf("progress: %d%%\n", renderProgress1);
                renderProgress0 = renderProgress1;
            }
            

        }
    }
    //printf("--------------------------------\n");
    for (unsigned int x = 0; x < WIDTH; x++) {
        for (unsigned int y = 0; y < HEIGHT; y++) {
            if ((float)render[y][x][0] > 255 || (float)render[y][x][1] > 255 || (float)render[y][x][2] > 255) {
                printf("somehow larger than 255\n");
            }
            if (render[y][x][0] < 0 || render[y][x][1] < 0 || render[y][x][2] < 0) {
                printf("somehow smaller than 0\n");
            }
        }
    }
  //display pixels according to buffer.
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      plot_pixel(x, y, render[y][x][0], render[y][x][1], render[y][x][2]);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

