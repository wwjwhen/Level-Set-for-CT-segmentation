# Level-Set-for-CT-segmentation
This is a level set method in functional form to segment CT images of main aortic

## Balance deformable and rigid
Set the negative parts of the level set function to zero if the corresponding pixel is below the threshold. 

## Make the *gradient* operator much faster with the help of *DGradient*

Find more details [here](https://cn.mathworks.com/matlabcentral/fileexchange/29887-dgradient?focused=5175102&tab=function)

The .mexw64 is generated on win64-vs2013

## Make the *del2* operator much faster by rewriting *del2* to *DEL2*

The .mexw64 is generated on win64-vs2013