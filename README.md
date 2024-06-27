# MLBAH88-48-88
%%%

DDDD     L          U     U     TTTTTTTTT
D   D    L          U     U         T
D    D   L          U     U         T
D    D   L          U     U         T
D   D    L          U     U         T
DDDD     LLLLLL      UUUUU          T        
                 
%%%

An 88-48-88 line MATLAB code for asymptotic homogenisation of spatially-varying multiscale configurations

## Getting Started

The code is documented in the paper: "An 88-48-88 line MATLAB code for asymptotic homogenisation of spatially-varying multiscale configurations, Chuang Ma, Shaoshuai Li, Dingchuan Xue, Yichao Zhu and Xu Guo, 2024".

The code was developed and tested using MATLAB, version R2021a, including MATLAB Deep Learning Toolbox.

The code package mainly contains three functions: ```Compute_DH.m```, ```NNtrain.m```, ```MLBAH.m```.

a) Compute_DH.m : an 88-line script to compute the homogenised properties of a given irregular-shaped building block;

b) NNtrain.m : a 48-line script to train a neural network expressing the homogenised properties;

c) MLBAH.m : an 88-line script to conduct online compliance computation.

Here are a few simple examples of calling this code package.

### Matlab example
```
% Calculate the effective moduli of the X-shaped unit cell with a volume fraction of 0.3.
load('...\Phi\phi_X.mat');
D_H=Compute_DH(phi_X,0.2958,[1 0;0 1],[1e-10 0.3297],[1e-10 0.3846]);

% Calculate the homogenised compliance of a 2D rectangular short beam (filled with gradient structure).
MLBAH(2,1,200,100,[ 1; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0 ],[ -0.06; -0.08; 0.07; -0.06; 0.19; 0.18 ]);
```
### Additional auxiliary functions
darkb2r.m : This is a color package function, copyright owned by Cunjie Zhang. Email: daisy19880411@126.com.
gradedstructure.m : This function is used for plotting spatially-varying multiscale configuration (SVMSC).
gradedstructure([ 1; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0 ],phi_X,[ -0.06; -0.08; 0.07; -0.06; 0.19; 0.18 ]);

## Help

Please send your comments or questions to: yichaozhu@dlut.edu.cn / chuangma@mail.dlut.edu.cn

## Authors

This Matlab code was written by Chuang Ma, Shaoshuai Li, Dingchuan Xue, Yichao Zhu and Xu Guo,
Dalian University of Technology,
116023 Dalian, China.                                                

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

The financial supports from National Key Research and Development Plan (2020YFB1709401) from the Ministry of Science and Technology of the People's Republic of China, and the National Natural Science Foundation of China (12172074).
