#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <map>
#include <random>

#include <glm/glm.hpp>

#include "utils.h"
#include "Object3D.h"
#include "../external_libraries/common_include/kdtree++/kdtree.hpp"

class Scene {


private:
    friend struct scene_traverser;
public:
    std::vector<Object3D*> object_list;
	std::vector<LightSource*> light_list;
    std::map<std::string, Material*> material_list;


	Scene(const char* file_path);
	~Scene();

    int getNumberOfTriangles();
	int getNumberOfObjects();
	int getNumberOfSpheres();

};

#endif 