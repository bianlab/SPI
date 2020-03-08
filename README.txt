--------------------------------------------------------------------------------------------------------------------------

Code demo for single-pixel imaging (SPI) with different reconstruction methods including
[1] differential ghost imaging (DGI)
[2] gradient descent (GD)
[3] conjugate gradient descent (CGD)
[4] Poisson maximum likelihood (Poisson)
[5] alternating projection (AP)
[6] sparse representation compressive sensing (Sparse)
[7] total variation compressive sensing (TV)


Public release v3.0 (Oct 24, 2017) 

------------------------------------------------------------------------------------------------------------------------------------
 Contents
------------------------------------------------------------------------------------------------------------------------------------
Note: Before running the code demo, please change Matlab current folder to ¡°../Code_SPI_pkg_3.0¡±.


***************** main function *********************

*) Demo.m                : This demo does the simulation of single-pixel imaging, with different reconstruction methods.

***************** reconstruction functions ********************

*) fun_SPI_R_DGI.m       : This function implements the singel-pixel imaging reconstruction using the differential ghost imaging (DGI) method.
W. Gong and S. Han. ¡®A method to improve the visibility of ghost images obtained by thermal light,¡¯ Phys. Lett. A, 374(8):1005-1008, 2010.
F. Ferri, D. Magatti, L.A. Lugiato, and A. Gatti. ¡®Differential ghost imaging,¡¯ Phys. Rev. Lett., 104(25):253603, 2010.

*) fun_SPI_R_GD.m        : This function implements the singel-pixel imaging reconstruction using the gradient descent method.
Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

*) fun_SPI_R_CGD.m       : This function implements the singel-pixel imaging reconstruction using the conjugate gradient descent method.
Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

*) fun_SPI_R_Poisson.m   : This function implements the singel-pixel imaging reconstruction using the Poisson maximum likelihood method.
Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

*) fun_SPI_R_AP.m        : This function implements the singel-pixel imaging reconstruction using the alternating projection method proposed in 
Kaikai Guo, Shaowei Jiang, and Guoan Zheng, ¡®Multilayer fluorescence imaging on a single-pixel detector,¡¯ Biomedical Optics Express, 7, 7, 2425 (2016).

*) fun_SPI_R_Sparse.m    : This function implements the singel-pixel imaging reconstruction using the compressive sensing method under sparse DCT (descrete cosine tranform) representation regulation.
Duarte, Marco F., et al. "Single-pixel imaging via compressive sampling." IEEE signal processing magazine 25.2 (2008): 83-91.
Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

*) fun_SPI_R_TV.m        : This function implements the singel-pixel imaging reconstruction using the compressive sensing method under total variation (TV) regulation.
Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.
Xuemei Hu, Jinli Suo, Tao Yue, Liheng Bian and Qionghai Dai, 'Patch-primitive driven compressive ghost imaging,'  Optics Express, 2015, vol. 23, no. 9, pp. 11092-11104.

------------------------------------------------------------------------------------------------------------------------------------
 Disclaimer
------------------------------------------------------------------------------------------------------------------------------------
Any unauthorized use of these routines for industrial or profit-oriented activities is expressively prohibited.

------------------------------------------------------------------------------------------------------------------------------------
 Feedback
------------------------------------------------------------------------------------------------------------------------------------
If this code offers help for your research, please cite our paper:
Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

If you have any comment, suggestion, or question, please feel free to contact Liheng Bian at lihengbian@gmail.com.
------------------------------------------------------------------------------------------------------------------------------------