#include "../include/Scenes.h"

#include "../external_libraries/common_include/pugixml.h"
#include "../include/xmlTraverser.h"

#include <iostream>
#include <string>
#include <random>


Scene::Scene (const char* file_path) {
	if (!file_path) {
		std::cout << "Error: no file path" << std::endl;
		exit (EXIT_FAILURE);
	}

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(file_path);
    std::cout << "xml result: " << result.description() << std::endl;

    if (!result) {
        std::cout << "Error: no xml result" << std::endl;
        exit (EXIT_FAILURE);
    }
		
	
	scene_traverser walker;
	walker.scene = this;

	doc.traverse(walker);
    std::cout << "Scene created!" << std::endl;
}

Scene::~Scene() {

    for (std::map<std::string, Material* >::iterator it = material_list.begin(); it != material_list.end(); it++) {
		delete it->second;
	}

    for (int i = 0; i < object_list.size(); ++i) {
        delete object_list[i];
    }

    for (int i = 0; i < light_list.size(); ++i) {
        delete light_list[i];
    }
	
}


int Scene::getNumberOfObjects() {
	return object_list.size();
}


int Scene::getNumberOfTriangles() {
	int num_triangles = 0;
	for (int i = 0; i < object_list.size(); ++i) {
		if (Mesh* m = dynamic_cast<Mesh*>(object_list[i]))
		   num_triangles += m->getNumberOfTriangles();
	}
	return num_triangles;
}



int Scene::getNumberOfSpheres()
{
	int num_spheres = 0;
	for (int i = 0; i < object_list.size(); ++i) {
		if (Sphere* s = dynamic_cast<Sphere*>(object_list[i]))
		    num_spheres++;
	}
	return num_spheres;
}

