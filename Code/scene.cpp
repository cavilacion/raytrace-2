#include "scene.h"

#include "hit.h"
#include "image.h"
#include "material.h"
#include "ray.h"

#include <cmath>
#include <limits>

using namespace std;

Color Scene::trace(Ray const &ray, const int depth)
{
		unsigned idx, obj_idx;
		
    // Find hit object and distance
    Hit min_hit(numeric_limits<double>::infinity(), Vector());
    ObjectPtr obj = nullptr;
    for (idx = 0; idx != objects.size(); ++idx)
    {
        Hit hit(objects[idx]->intersect(ray));
        if (hit.t < min_hit.t)
        {
						obj_idx = idx;
            min_hit = hit;
            obj = objects[idx];
        }
    }

    // No hit? Return background color.
    if (!obj) return Color(0.0, 0.0, 0.0);
		
    Material material = obj->material;             //the hit objects material
    Point hit = ray.at(min_hit.t);                 //the hit point
    Vector N = min_hit.N;                          //the normal at hit point
    Vector V = -ray.D;                             //the view vector
    Vector I;																			 //the total color
    Color color = material.color;
		
    
		for(auto lightptr : lights) {
			
			Triple Ia = color*material.ka;  
			Vector lightposition = lightptr->position;
			Color lightcolor = lightptr->color;
			Vector lightdir = (hit-lightposition).normalized();
			
			// check if light source hits other object first
			int in_shadow = 0; // assume object is not in shadow
			if (shadows) for (idx = 0; idx < objects.size(); ++idx) {
				if (idx != obj_idx) {
					Ray lightray(lightposition, lightdir);
					Hit this_hit(objects[obj_idx]->intersect(lightray));
					Hit other_hit(objects[idx]->intersect(lightray));
					if (other_hit.t < this_hit.t) {
							in_shadow = 1; 
							break;
					}
				}
			}
			
			if (!in_shadow) {
				Vector L = (lightposition-hit).normalized();
				Vector R = (2.0*N.dot(L)*N) - L;			// Reflection fector 
												
				Triple Id = color*lightcolor*material.kd*max(0.0,L.dot(N));
				Triple Is = pow(max(0.0,V.dot(R)),material.n)*lightcolor*material.ks;
				
				I += Id + Is;
			}
			
			I += Ia;
		}   
		
		// continuation of view ray (specular reflection) with recursive call
		if (depth < num_reflections) {
			Vector reflect = (2.0*N.dot(V)*N - V).normalized();
			Ray from(hit+0.05*reflect, reflect);
			Color light_reflection = material.ks*trace (from, depth+1);
			I += light_reflection;
		}
		return I ;
}

void Scene::render(Image &img)
{
    unsigned w = img.width();
    unsigned h = img.height();
    for (unsigned y = 0; y < h; ++y)
    {
        for (unsigned x = 0; x < w; ++x)
        {
						int num_pixels = supersampling*supersampling;
						float start = 0.5/supersampling;
						float delta = 1.0/supersampling;
						Color *colors = (Color*)malloc (num_pixels*sizeof(Color));
						for (int i = 0; i < supersampling; i++) {
								for (int j = 0; j < supersampling; j++) {
										Point pixel(x + (start + i*delta), h - 1 - y + (start + j*delta), 0);
										Ray ray(eye, (pixel - eye).normalized());
										Color col = trace(ray, 0);
										col.clamp();
										colors[i*supersampling + j] = col;
								}
						}
						Color col;
						float weight = 1.0/num_pixels;
						for (int i = 0; i < num_pixels; i++) {
								col += weight*colors[i];
						}
						img(x, y) = col;
        }
    }
}

// --- Misc functions ----------------------------------------------------------

void Scene::addObject(ObjectPtr obj)
{
    objects.push_back(obj);
}

void Scene::addLight(Light const &light)
{
    lights.push_back(LightPtr(new Light(light)));
}

void Scene::setEye(Triple const &position)
{
    eye = position;
}

void Scene::configure(bool show_shadows, int num_ref, int supsam) {
	shadows = show_shadows;
	num_reflections = num_ref;
	supersampling = supsam;
}

unsigned Scene::getNumObject()
{
    return objects.size();
}

unsigned Scene::getNumLights()
{
    return lights.size();
}
