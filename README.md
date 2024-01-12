# Xplot
Xplot is a small Python command-line program to render images of atomic structures using [ASE](https://wiki.fysik.dtu.dk/ase/index.html) and [POV-Ray](http://www.povray.org/).

## Installation
First, you need to download the program by cloning this repository into your local machine:  
```sh
git clone https://github.com/Enry99/xplot.git
```
Once the download is completed, go inside the downloaded folder and run: 

```sh
bash install.sh
```

to add the executable in the `PATH` variable. With this operation, the you will be able to launch the program from any folder. The next step is to install the required dependencies stored in the 'requirements.txt' file. This can be performed with the command:

```sh
pip install -r requirements.txt
```

In order to generate the images and the movies, you need the *povray* and [*ffmpeg*](https://ffmpeg.org/) packages, which can be installed on Ubuntu with:

```sh
sudo apt install povray ffmpeg
```

## Usage

Basic usage for generating an image of the (final) configuration:
```sh
xplot filename.xyz
```

To generate the frames for a trajectory you can use:
```sh
xplot filename.xyz -i :
```

To plot only a frame every 10 frames you can use:
```sh
xplot filename.xyz -i ::10
```

To generate also a movie with the frames you can use:
```sh
xplot filename.xyz -i ::10 -m
```

A list with all the options can be found with 
```sh
xplot -h
```

**Important note**: the image generation with *povray* only works in Linux. For Windows, you can still use this program without povray by adding the option `-np`. However, the image quality is worse and bonds are not drawn.
