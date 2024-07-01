# JLSD - Julia SerDes
This repository hosts Pluto notebooks and relevant source codes for a personal project - using Julia to build fast and lightweight SerDes models and simulation framework

The goals for this project (still in its infant stage) are
- Document and share my personal Julia journey so far
- Using SerDes simulation as an example to demonstrate the pros/cons of Julia
- Begin a open-source and expandable SerDes simulation framework that encourages academia and industry adoption to evaluate more sophositicated architectures and algorithms
- Expose my fellow analog/mixed-signal designers to Julia (because not everyone can get a MATLAB license). It's much easier to design circuits when one can play with the specifications instead of taking it from others at face value

## Simulation framework
The modeling framework is based on custom structs and block simulations. The code is not heavily documentated yet, but should be self-explanatory. Go through the Pluto notebooks to understand the key concepts in the models.

## Notebooks
In the Pluto Notebooks directory, you will find the .jl files for the notebooks to be viewed and played with on your local machine. .html and .pdf files are also included in the directory.

## Standalone widget
Currently, there is a demo widget powered by Makie (shown below. Simply run Main_UI.jl and start playing with a basic SerDes model's parameters. The model consists of a relatively detailed transmitter, a low-loss channel, and a simple receiver w/ baud-rate CDR. Note that the widget might be continuously updated to include more (common) features. Use this as an example to extend to your own models. 
![widget_ui](https://github.com/kevjzheng/JLSD/blob/main/img/widget_ui.png)
