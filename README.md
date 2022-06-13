# CG-project

### **Title**

Monte Carlo Ray Tracer

### **Summary**

A program for rendering virtual 3D scenes with the optimized Monte Carlo ray tracing.

### **Team Members**

Ruofan Xu, Xinyi Shi, Huaijin Wu

### **Problem Description**

When rendering images we often strive for photo realism. In order to solve this, we need to solve the global illumination problem, which means simulating how light behaves in the real world. In this project, we'll use Monte Carlo ray tracing algorithm to simulate many light phenomena, such as color bleeding, soft shadows, reflection and refraction.

The challenge lies in building a rendering system from scratch, from which we can learn a big picture of rendering and get more familiar with global illumination. Through implementing different techniques, we can learn how to improve the quality of rendered images step by step.

A key tool for getting many effects accurately in a ray tracer is Monte-Carlo integration. (We will solve this by learning) So we will learn basic monte carlo methods for ray tracing, then read other related papers and theories to implement further techniques to get better effects or shorten the running time.

### **Goals and Deliverables**

1）Implement a cpu(maybe GPU) renderer using Monte-Carlo ray tracing algorithm

2）Use KDTree to accelerate(different from the BVH we have implemented in scotty3d)

3）Implement parallelization using openMP

4）Implement XML parser to load XML files

5）Implement more graphic techniques (one or more ):

1. photon mapping
2. Make different wavelengths refract a different amount and increase the number of wavelengths.(this can achieve rainbow effects from refraction)
3. micro facet model

Our baseline plan: implement a renderer using Monte Carlo ray tracing algorithm to at least get the image as well as our scotty3d. What we expect is to have a higher quality and render faster. We also hope our system can achieve some fancy contents, including rainbow effects from refraction and micro facet. We evaluate our results by the scene we finally get.

### **Schedule**

week 1: Read related websites and papers and learn about detailed theories of Monte-Carlo ray tracing algorithm

week 2: Build the basic framework of code 

week 3: Realize the framework of code

week 4: Debug and test the code/see how it works

week 5: Discuss the ways to optimize the results and implement some optimizations

week 6: Draw a conclusion of our work. Make the video, write our paper and prepare for the presentation.

### **Building**

```
git clone https://github.com/wuhuaijin/CG-project
cd CG-project
mkdir build
cd build
cmake ..
make
```

This will build files in the folder `build` of the cloned repository, which can be used to build the program. It mail fail if version of g++ or gcc dismatch, please check 4th and 5th lines of `CMakeLists.txt`.

### **Usage**

```
./global_illumination ../data/scenes/xxxx.xml
```

After rendering is finished, you can find the generated xml file in folder `build`.

### **ResourcesList**

- [Physically Based Rendering](http://www.pbr-book.org/) - Matt Pharr, Wenzel Jakob and Greg Humphreys

- [Global Illumination using Photon Maps](http://graphics.stanford.edu/~henrik/papers/ewr7/ewr7.html) - Henrik Wann Jensen
- https://www.pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models
- http://www.cs.cornell.edu/courses/cs4620/2013fa/lectures/22mcrt.pdf
