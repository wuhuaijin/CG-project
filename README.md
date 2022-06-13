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
