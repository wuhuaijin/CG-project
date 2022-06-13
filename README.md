# CG-project

### **Title**

Monte Carlo Ray Tracer

### **Summary**

A program for rendering virtual 3D scenes with the optimized Monte Carlo ray tracing.

### **Team Members**

Ruofan Xu, Xinyi Shi, Huaijin Wu

### **Building**

Parameters in `main.cpp` can be modified in order to set resolution and subsampling . Then build the excutable:

```
git clone https://github.com/wuhuaijin/CG-project
cd CG-project
mkdir build
cd build
cmake ..
make
```

This will build files in the folder `build` of the cloned repository, which can be used to build the program. Cmake mail fail if version of g++ or gcc dismatches `CMakeLists.txt`, please check 4th and 5th lines of `CMakeLists.txt`.

### **Usage**

Use the excutable you generated in folder `build`, and choose the scene you want to render:

```
./global_illumination ../data/scenes/xxxx.xml
```

After rendering is finished, you can find the generated xml file in the folder `build`.

OR you can simply use the excutable we provide(in the folder `results`). 

### **ResourcesList**

- [Physically Based Rendering](http://www.pbr-book.org/) - Matt Pharr, Wenzel Jakob and Greg Humphreys

- [Global Illumination using Photon Maps](http://graphics.stanford.edu/~henrik/papers/ewr7/ewr7.html) - Henrik Wann Jensen
- https://www.pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models
- http://www.cs.cornell.edu/courses/cs4620/2013fa/lectures/22mcrt.pdf
