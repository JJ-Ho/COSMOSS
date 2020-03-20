What is COSMOSS?
===
:::info
**COSMOSS** is a **C**oupled **OS**cillator **MO**del **S**pectral **S**imulator
:::


<figure style="float:right;">
    <img src="https://i.imgur.com/HOQWa4m.png" align="right" width="160" height="160">
</figure>

I originally wrote COSMOSS for my own research projects. Using the coupled-oscillators framwork, I simulate and extract useful information from spectrum in my [previous works](https://scholar.google.com/citations?user=O8Y6wd8AAAAJ&hl=en&authuser=1), mostly on 2D IR/ 2D SFG. As the COSMOSS project moving forward, I couldn't help but wonder if I can make it more general (as its name implied), i.e. works for different molecules, with a user-freindly GUI, support other spectral types, and most important of all, it is **open source**! So here is it, hope you enjoy it.   

## What can COSMOSS do?
COSMOSS can help you generate model structures and simulated their vibrational spectrum, including:

- Fourier Transform Infrared spectrum [(FTIR)](https://en.wikipedia.org/wiki/Fourier-transform_infrared_spectroscopy)
- Sum-Frequency Generation spectrum [(SFG)](https://en.wikipedia.org/wiki/Sum_frequency_generation_spectroscopy)
- Two-Dimensional Infrared spectrum [(2D IR)](https://en.wikipedia.org/wiki/Two-dimensional_infrared_spectroscopy)
- Two-Dimensional Sum-Frequency Generation spectrum [(2D SFG)](https://zanni.chem.wisc.edu/2d-sfg-spectroscopy/)

Besides simple spectral simulations, you can also use it in a creative way. For example, I used COSMOSS to study the interaction between peptide aggregates and self-assembly monolayers and found 2D SFG is very informative! [^Ref_Paper_1] 
![](https://i.imgur.com/5M4loNP.gif)

[^Ref_Paper_1]: [J. Phys. Chem. A 2018, 122, 5, 1270-1282](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.7b11934)

The [coding structure](/Vwfo8DgmQL2Ya_ryUYd06w) is very general! As long as the 'modes' in question can be model as  coupled oscillators, you can generate their corresponding spectrum as well. For example, if you are interested in simulating the UV-Vis spectra of a model system, you can repurpose the FTIR core with little modification. 


Documentation
---
![](https://i.imgur.com/pDJrSvz.jpg)

All the notes and tutorails are hosted on the HackMD. Learn more [in the documentation](https://hackmd.io/@jjh/COSMOSS) This documentation will help you understand the coding structure and theoretical backgroung so you can make the best of COSMOSS.




Installation
---
[![View COSMOSS on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/64433-cosmoss)

Download the source code from either [Matlab file-exchange](https://www.mathworks.com/matlabcentral/fileexchange/64433-cosmoss) or [GitHub](https://github.com/JJ-Ho/COSMOSS).[^Download] Once download, double click on the 'COSMOSS.mlapp' will bring up the main GUI! [^mlapp]

Collaboration
---
One of the main reason why COSMOSS is open source is for collaboration! Any contribution is welcomed. The are two branches for potential collaboration:

- Documentation
- Source code

This website is editable if you regist an account and sign-in.  Lean more on how to edit the page on HackMD here: 
- [HeckMD help page (English)](https://hackmd.io/c/tutorials/%2Fs%2Ftutorials). 
- [HeckMD help page (Chinese)](https://hackmd.io/c/tutorials-tw/%2Fs%2Ftutorials-tw). 

Feel free to write/edit about anything relevent. You are also welcome to send me the link of your personal documentation though my email: 

- jhoatgit [at] gmail.com

or leave a comment on this site. 

For the source code, please use the [GitHub](https://github.com/JJ-Ho/COSMOSS).


*[COSMOSS]: Coupled OScillator MOdel Spectral Simulator
*[FTIR]: Fourier Transform Infrared spectrum
*[UV-Vis]: UV-Visible spectroscopy

[^Download]:Since v2.3.1, the [Matlab file-exchange](https://www.mathworks.com/matlabcentral/fileexchange/64433-cosmoss) download page is integrated with the master branch of the [GitHub](https://github.com/JJ-Ho/COSMOSS) repository. No more fooling around!

[^mlapp]: You can modify the GUI elements and the callbacks with the [AppDesinger](https://www.mathworks.com/videos/app-designer-overview-1510748719083.html) in Matlab. See [How to edit COSMOSS GUI](/Vy7kwQ4WTHGCOLoZqNPhnw?sync=&type=) for more information.

###### tags: `COSMOSS`



