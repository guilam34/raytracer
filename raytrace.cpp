//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

int g_width;
int g_height;

struct Ray
{
	vec4 origin;
	vec4 dir;
};

//Sphere structure to hold sphere data
struct Sphere
{
	string name;
	vec4 pos;
	vec3 scl;
	vec3 rgb;
	float k_amb;
	float k_diff;
	float k_spec;
	float k_refl;
	float spec_comp;
};

//Light structure to hold light data
struct Light
{
	string name;
	vec4 pos;
	vec3 intensity;
};


//Sphere vector to hold all spheres in input file
vector<Sphere> list_spheres;

//Light vector to hold all lights in input file
vector<Light> list_lights;

vector<vec4> g_colors;
vec4 back_rgb;
vec4 ambient_rgb;
float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

string output_file;

// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
	stringstream ss(s1 + " " + s2 + " " + s3);
	vec4 result;
	ss >> result.x >> result.y >> result.z;
	result.w = 1.0f;
	return result;
}

float toFloat(const string& s)
{
	stringstream ss(s);
	float f;
	ss >> f;
	return f;
}

void parseLine(const vector<string>& vs)
{
	/*Parsing input file line by line; each line is a vector of strings with the 0th index
	denoting the field and the ensuing elements being the values associated for the field*/
	if (vs[0] == "RES")
	{
		g_width = (int)toFloat(vs[1]);
		g_height = (int)toFloat(vs[2]);
		g_colors.resize(g_width * g_height);
	}
	else if (vs[0] == "NEAR")
	{
		g_near = -1.*(int)toFloat(vs[1]);
	}
	else if (vs[0] == "LEFT")
	{
		g_left = (int)toFloat(vs[1]);
	}
	else if (vs[0] == "RIGHT")
	{
		g_right = (int)toFloat(vs[1]);
	}
	else if (vs[0] == "BOTTOM")
	{
		g_bottom = (int)toFloat(vs[1]);
	}
	else if (vs[0] == "TOP")
	{
		g_top = (int)toFloat(vs[1]);
	}
	else if (vs[0] == "SPHERE")
	{
		/* We stored scale, position, and rgb of the sphere in vectors to better organize
		our sphere structure. Each new sphere read is pushed onto a vector of spheres*/
		Sphere newSphere;
		newSphere.name = (string)vs[1];
		newSphere.pos = vec4(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), 1.0f);
		newSphere.scl = vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
		newSphere.rgb = vec3(toFloat(vs[8]), toFloat(vs[9]), toFloat(vs[10]));
		newSphere.k_amb = toFloat(vs[11]);
		newSphere.k_diff = toFloat(vs[12]);
		newSphere.k_spec = toFloat(vs[13]);
		newSphere.k_refl = toFloat(vs[14]);
		newSphere.spec_comp = toFloat(vs[15]);
		list_spheres.push_back(newSphere);
	}
	else if (vs[0] == "BACK")
	{
		back_rgb = vec4(toFloat(vs[1]), toFloat(vs[2]), toFloat(vs[3]), 1.0f);
	}
	else if (vs[0] == "OUTPUT")
	{
		output_file = vs[1];
	}
	else if (vs[0] == "LIGHT")
	{
		Light newLight;
		newLight.name = vs[1];
		newLight.pos = vec4(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), 1.0f);
		newLight.intensity = vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
		list_lights.push_back(newLight);
	}
	else if (vs[0] == "AMBIENT")
	{
		ambient_rgb = vec3(toFloat(vs[1]), toFloat(vs[2]), toFloat(vs[3]));
	}
}

void loadFile(const char* filename)
{
	ifstream is(filename);
	if (is.fail())
	{
		cout << "Could not open file " << filename << endl;
		exit(1);
	}
	string s;
	vector<string> vs;
	while (!is.eof())
	{
		vs.clear();
		getline(is, s);
		istringstream iss(s);
		while (!iss.eof())
		{
			string sub;
			iss >> sub;
			vs.push_back(sub);
		}
		parseLine(vs);
	}
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
	int iy2 = g_height - iy - 1; // Invert iy coordinate.
	g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine
vec4 trace(const Ray& ray, int trace_mode, float light_dist, int& num_recurrences);
vec4 pixel_lighting(vec4 closest_intpoint, vec4 origin, mat4 inv_sphereMat, Sphere curSphere, vec4 unit_normal, const Ray& ray, int num_recurrences)
{
	//Return value
	vec4 lighting(1, 1, 1, 1);

	//Light reflection ray declaration
	vec4 refl_vec;

	//Variables to store shadow ray calculations
	Light curLight;
	vec4 light_vec;
	float light_veclen;
	vec4 norm_light;

	//Calculation of normalized vector v for specular effect
	vec4 v = origin - closest_intpoint;
	float v_len = sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
	vec4 norm_v = vec4(1.*v.x / v_len, 1.*v.y / v_len, 1.*v.z / v_len, 0.0f);

	//Normal calculation for the closest intersection point for a ray by transposing normal vector on unit sphere	
	vec4 trans_normal = transpose(inv_sphereMat)*unit_normal;
	float trans_len = sqrt(1.*pow(trans_normal.x, 2) + 1.*pow(trans_normal.y, 2) + 1.*pow(trans_normal.z, 2));
	trans_normal = vec4(1.*trans_normal.x / trans_len, 1.*trans_normal.y / trans_len, 1.*trans_normal.z / trans_len, 0.0f);

	//Reflection effect calculations
	vec4 refl_effect(0, 0, 0, 0);
	if (num_recurrences <= 3 && curSphere.k_refl != 0.0f)
	{
		vec4 refl_v = -2 * (trans_normal.x*ray.dir.x + trans_normal.y*ray.dir.y + trans_normal.z*ray.dir.z)*trans_normal + ray.dir;
		refl_v = normalize(refl_v);
		Ray refl_ray;
		refl_ray.origin = closest_intpoint;
		refl_ray.dir = refl_v;
		refl_effect = trace(refl_ray, 0, 0, num_recurrences);
		if (refl_effect.x == back_rgb.x && refl_effect.y == back_rgb.y && refl_effect.z == back_rgb.z)
		{
			refl_effect = vec4(0.0, 0.0, 0.0, 1.0f);
		}
	}

	//Vectors to store ambient, diffuse, and specular effects of all light sources
	vec4 diffuse_effect(0, 0, 0, 0);
	vec4 specular_effect(0, 0, 0, 0);
	vec4 ambient_effect = vec4(ambient_rgb.x*curSphere.k_amb*curSphere.rgb.x,
		ambient_rgb.y*curSphere.k_amb*curSphere.rgb.y,
		ambient_rgb.z*curSphere.k_amb*curSphere.rgb.z, 1.0f);

	//Iterate through every light in light vector
	vector<Light>::iterator light_it = list_lights.begin();
	for (light_it; light_it != list_lights.end(); light_it++)
	{
		//Calculate normalized shadow ray
		curLight = (*light_it);
		light_vec = curLight.pos - closest_intpoint;
		light_veclen = sqrt(pow(light_vec.x, 2) + pow(light_vec.y, 2) + pow(light_vec.z, 2));
		norm_light = vec4(1.*light_vec.x / light_veclen, 1.*light_vec.y / light_veclen, 1.*light_vec.z / light_veclen, 0.0f);

		//Initialize a ray to test light-sphere intersections
		Ray intersection_ray;
		intersection_ray.dir = curLight.pos - closest_intpoint;
		intersection_ray.dir = normalize(intersection_ray.dir);
		intersection_ray.origin = closest_intpoint;

		//Check if light ray intersects any spheres
		int light_free = 0;
		int rec_fake = 0;
		vec4 trace_retval = trace(intersection_ray, 1, length(curLight.pos - closest_intpoint), rec_fake);
		if (trace_retval.x == back_rgb.x && trace_retval.y == back_rgb.y && trace_retval.z == back_rgb.z)
		{
			light_free = 1;
		}

		//Diffuse effect calculation for ONE light source
		float normTimesLight = (1.*trans_normal.x*norm_light.x) + (1.*trans_normal.y*norm_light.y) + (1.*trans_normal.z*norm_light.z);
		if (normTimesLight >= 0 && light_free == 1)
		{
			diffuse_effect = vec4(diffuse_effect.x + (curLight.intensity.x*curSphere.k_diff*normTimesLight*curSphere.rgb.x),
				diffuse_effect.y + (curLight.intensity.y*curSphere.k_diff*normTimesLight*curSphere.rgb.y),
				diffuse_effect.z + (curLight.intensity.z*curSphere.k_diff*normTimesLight*curSphere.rgb.z), 1.0f);
		}
		else{
			continue;
		}
		//Calculate light reflection ray
		refl_vec = 2 * trans_normal*(trans_normal.x*norm_light.x + trans_normal.y*norm_light.y + trans_normal.z*norm_light.z) - norm_light;

		//Specular effect calculation for ONE light source
		float rTimesV = refl_vec.x*norm_v.x + refl_vec.y*norm_v.y + refl_vec.z*norm_v.z;
		if (rTimesV >= 0 && light_free == 1)
		{
			specular_effect = vec4(specular_effect.x + (curLight.intensity.x*curSphere.k_spec*pow(rTimesV, curSphere.spec_comp)),
				specular_effect.y + (curLight.intensity.y*curSphere.k_spec*pow(rTimesV, curSphere.spec_comp)),
				specular_effect.z + (curLight.intensity.z*curSphere.k_spec*pow(rTimesV, curSphere.spec_comp)), 1.0f);
		}
		else {
			continue;
		}
	}

	//Sum up all lighting components for a single pixel
	lighting = vec4(diffuse_effect.x + specular_effect.x + ambient_effect.x + curSphere.k_refl*refl_effect.x,
		diffuse_effect.y + specular_effect.y + ambient_effect.y + curSphere.k_refl*refl_effect.y,
		diffuse_effect.z + specular_effect.z + ambient_effect.z + curSphere.k_refl*refl_effect.z,
		1.0f);

	return lighting;
}

// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int trace_mode, float light_dist, int& num_recurrences)
{
	//Reflection recurrences break condition
	if (trace_mode == 0 && num_recurrences == 3)
		return vec4(0, 0, 0, 0.0f);

	//Increment recurrence counter only if we're not checking light intersection
	if (trace_mode == 0)
	{
		num_recurrences = num_recurrences + 1;
	}
	//Declare vec4 variables to hold transformed ray 
	mat4 inv_sphereMat;
	vec4 trans_dir;
	vec4 trans_origin;
	vec4 pixel_rgb;
	vec4 lighting;
	vec4 int_point;
	vec4 closest_intPoint;
	vec4 unit_normal;
	bool found_inter = 1;
	float closest_z = 5000;
	float a;
	float b;
	float c;
	float trans_dirlen;
	float trans_originlen;
	float time1;
	float time2;

	//Iterate through vector of spheres finding
	vector<Sphere>::iterator sphere_it = list_spheres.begin();
	for (sphere_it; sphere_it != list_spheres.end(); sphere_it++)
	{
		//Generate sphere transformation matrix
		mat4 sphereMat(vec4(1.*(*sphere_it).scl.x, 0.0f, 0.0f, 1.*(*sphere_it).pos.x),
			vec4(0.0f, 1.*(*sphere_it).scl.y, 0.0f, 1.*(*sphere_it).pos.y),
			vec4(0.0f, 0.0f, 1.*(*sphere_it).scl.z, 1.* (*sphere_it).pos.z),
			vec4(0.0f, 0.0f, 0.0f, 1.0f));

		//Invert sphere transformation matrix so we can apply it to the ray instead
		InvertMatrix(sphereMat, inv_sphereMat);
		trans_dir = inv_sphereMat*ray.dir;
		trans_origin = inv_sphereMat*ray.origin;

		/*Computations of the coefficients for the quadratic equation used to find
		sphere-ray intersection times*/
		trans_dirlen = 1.*sqrt(1.*pow(1.*trans_dir.x, 2) + 1.*pow(1.*trans_dir.y, 2) + 1.*pow(1.*trans_dir.z, 2));
		trans_originlen = 1.*sqrt(1.*pow(1.*trans_origin.x, 2) + 1.*pow(1.*trans_origin.y, 2) + 1.*pow(1.*trans_origin.z, 2));
		a = 1.*pow(1.*trans_dirlen, 2);
		b = (1.*trans_dir.x*trans_origin.x) + (1.*trans_dir.y*trans_origin.y) + (1.*trans_dir.z*trans_origin.z);
		c = 1.*pow(1.*trans_originlen, 2) - 1;
		float quad_result = 1.*pow(1.*b, 2) - (1.*a*c);

		/*Use results of determinant computation of sphere-ray intersection equation
		to find sphere-ray intersection points*/
		if (quad_result <0)
		{
			continue;
		}

		//If there is only one intersection point then this case is fulfilled
		else if (quad_result == 0)
		{
			time1 = (-b*1.) / (1.*a);

			//For camera-sphere intersections: Check if intersection point is z-closest found and compute pixel lighting
			if (trace_mode == 0)
			{
				int_point = ray.origin + (time1*ray.dir);
				if (abs(int_point.z - ray.origin.z) < closest_z &&(time1 > 1||(num_recurrences>1&&time1>0.0001)))
				{
					closest_z = int_point.z;
					closest_intPoint = int_point;
					unit_normal = trans_origin + time1*trans_dir;
					lighting = pixel_lighting(closest_intPoint, ray.origin, inv_sphereMat, (*sphere_it), vec4(unit_normal.x, unit_normal.y, unit_normal.z, 0.0f), ray, num_recurrences);
					pixel_rgb = lighting;
					found_inter = 0;
				}
			}

			//For point-light intersections: Notify caller that intersection was found for light ray
			else if (trace_mode == 1)
			{
				vec4 light_intvec = time1*ray.dir;
				if (length(light_intvec) < light_dist && time1>0.001)
				{
					return vec4(0, 0, 0, 0.0f);
				}
			}
		}

		//If there are two intersection points then this case is fulfilled
		else if (quad_result > 0)
		{
			//Same process as for quad_result==0
			time1 = ((1.*-b) / (1.*a)) + (sqrt(1.*pow(1.*b, 2) - (1.*a*c)) / (1.*a));
			if (trace_mode == 0)
			{
				int_point = ray.origin + (time1*ray.dir);
				if (abs(int_point.z - ray.origin.z) < closest_z && (time1 > 1 || (num_recurrences>1 && time1>0.0001)))
				{
					closest_z = abs(int_point.z - ray.origin.z);
					closest_intPoint = int_point;
					unit_normal = trans_origin + time1*trans_dir;
					lighting = pixel_lighting(closest_intPoint, ray.origin, inv_sphereMat, (*sphere_it), vec4(unit_normal.x, unit_normal.y, unit_normal.z, 0.0f), ray, num_recurrences);
					pixel_rgb = lighting;
					found_inter = 0;
				}
			}
			else if (trace_mode == 1)
			{
				vec4 light_intvec = time1*ray.dir;
				if (length(light_intvec) < light_dist && time1>0.001)
				{
					return vec4(0, 0, 0, 0.0f);
				}
			}

			//Repeat the process for the second intersection
			time2 = ((1.*-b) / (1.*a)) - (sqrt(1.*pow(1.*b, 2) - (1.*a*c)) / (1.*a));
			if (trace_mode == 0)
			{
				int_point = ray.origin + (time2*ray.dir);
				if (abs(int_point.z - ray.origin.z) < closest_z && (time2 > 1 || (num_recurrences>1 && time2>0.0001)))
				{
					closest_z = abs(int_point.z - ray.origin.z);
					closest_intPoint = int_point;
					unit_normal = trans_origin + time2*trans_dir;
					lighting = pixel_lighting(closest_intPoint, ray.origin, inv_sphereMat, (*sphere_it), vec4(unit_normal.x, unit_normal.y, unit_normal.z, 0.0f), ray, num_recurrences);
					pixel_rgb = lighting;
					found_inter = 0;
				}
			}
			else if (trace_mode == 1)
			{
				vec4 light_intvec = time2*ray.dir;
				if (length(light_intvec) < light_dist && time2>0.001)
				{
					return vec4(0, 0, 0, 0.0f);
				}
			}
		}
	}
	//If there was no intersection found return the background color
	if (found_inter == 1)
		return back_rgb;
	else
	{
		return pixel_rgb;
	}
}

vec4 getDir(int ix, int iy)
{
	//Get the x,y camera coordinates of the pixels first with camera coordinates pointing to lower left of a pixel
	float real_x = g_left + (1.*ix / (1.*g_width))*(g_right - g_left);
	float real_y = g_bottom + (1.*iy / (1.*g_height))*(g_top - g_bottom);
	return vec4(real_x, real_y, g_near, 0.0f);
}

void renderPixel(int ix, int iy)
{
	Ray ray;
	ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	ray.dir = getDir(ix, iy);
	int num_recurrences = 0;
	vec4 color = trace(ray, 0, -1, num_recurrences);
	setColor(ix, iy, color);
}

void render()
{
	for (int iy = 0; iy < g_height; iy++)
		for (int ix = 0; ix < g_width; ix++)
			renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels)
{
	FILE *fp;
	const int maxVal = 255;

	printf("Saving image %s: %d x %d\n", fname, Width, Height);
	fp = fopen(fname, "wb");
	if (!fp) {
		printf("Unable to open file '%s'\n", fname);
		return;
	}
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);

	for (int j = 0; j < Height; j++) {
		fwrite(&pixels[j*Width * 3], 3, Width, fp);
	}

	fclose(fp);
}

void saveFile()
{
	// Convert color components from floats to unsigned chars.
	// TODO: clamp values if out of range.
	unsigned char* buf = new unsigned char[g_width * g_height * 3];
	for (int y = 0; y < g_height; y++)
		for (int x = 0; x < g_width; x++)
			for (int i = 0; i < 3; i++)
			{
				if (g_colors[y*g_width + x][i] > 1.0)
					g_colors[y*g_width + x][i] = 1.0f;
				if (g_colors[y*g_width + x][i] < 0.0)
					g_colors[y*g_width + x][i] = 0.0f;
				buf[y*g_width * 3 + x * 3 + i] = (unsigned char)(((float*)g_colors[y*g_width + x])[i] * 255.9f);
			}

	//Get char pointer of output file name to pass into savePPM
	char* cstr = new char[output_file.length() + 1];
	strcpy(cstr, output_file.c_str());
	savePPM(g_width, g_height, cstr, buf);
	delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
	/*if (argc < 2)
	{
		cout << "Usage: template-rt <input_file.txt>" << endl;
		exit(1);
	}*/
	loadFile("testIllum.txt"); //argv[1]
	render();
	saveFile();
	return 0;
}

