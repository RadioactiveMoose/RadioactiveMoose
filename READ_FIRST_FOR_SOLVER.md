## Hello dear reader !

This work is still ongoing, thus here are some bullet points :

## Can I use this module for my research ?

That is what opensource is about ! So as long as I am credited, please do !
Crediting me allows me easier request for developping the solver and generally asking for help !

Moreover, do not hesitate to feedback on the solver, whether it's a chat, thank or bug report ;) 
Feel free to contact me at guillem.khairy@uni.lu or khairy.guillem@gmail.com

## What is the purpose of this module ?

The solver aims to study boiling and associated heat transfer phenomena in microgravity

## Which version of OpenFOAM ?

The OpenFOAM version used for this solver is OpenFOAM 12 (Foundation), porting to OpenFOAM 13 will happen soon (October-November 2025). 
Normally it is a "copy-paste" task but having not watched closely yet I will not warranty that it is.

## What boiling models are availables now ?

Right now only the Lee boiling model is featured and it is not benchmarked. 
When done, I will publish a table of the tested fluids and associated parameters. 

## Will other boiling model be available ?

Following the structure of fvModels/VoFCavitation, it is possible to add more boiling models. 
I intend to do so in the next months with the Schrage models and other relatively simple models.
RPI might come but is not a priority as for now. 

## Other addition planned to the module ?

I intend to also implement the Shear Lift Force as fvModels before beginning 2026.

## References :

1) R. Sugrue, J. Buongiorno, "A modified force-balance model for prediction of bubble departure diameter in subcooled flow boiling",Nuclear Engineering and Design, Volume 305,2016,Pages 717-722, ISSN 0029-5493, https://doi.org/10.1016/j.nucengdes.2016.04.017. (https://www.sciencedirect.com/science/article/pii/S0029549316300450)
2) Jeongmin Lee, Lucas E. O'Neill, Issam Mudawar, 3-D computational investigation and experimental validation of effect of shear-lift on two-phase flow and heat transfer characteristics of highly subcooled flow boiling in vertical upflow, International Journal of Heat and Mass Transfer, Volume 150, 2020, 119291, ISSN 0017-9310, https://doi.org/10.1016/j.ijheatmasstransfer.2019.119291.(https://www.sciencedirect.com/science/article/pii/S0017931019343741)
3) Jeongmin Lee, Issam Mudawar, Mohammad M. Hasan, Henry K. Nahra, Jeffrey R. Mackey, Experimental and computational investigation of flow boiling in microgravity, International Journal of Heat and Mass Transfer, Volume 183, Part C, 2022, 122237, ISSN 0017-9310, https://doi.org/10.1016/j.ijheatmasstransfer.2021.122237. (https://www.sciencedirect.com/science/article/pii/S0017931021013363)

## Acknowledgements :

Many thanks to the OpenFOAM community making opensource and free access to knowledge and science easier by the day. 
