# render_atoms
Small Python command-line program to render images of atomic structures using ASE and POVray.

To add the script to PATH, run "  bash install.sh  "

In order to generate the images and the movies, you need the 'povray' and 'ffmpeg' packages, which can be installed on Ubuntu with: \
" sudo apt install povray ffmpeg  "

A basic usage for generating an image of the (final) configuration is:
"xplot filename.xyz"

To generate the frames for a trajectory you can use:
"xplot filename.xyz -i :"

To plot only a frame every 10 frames you can use:
"xplot filename.xyz -i ::10"

To generate also a movie with the frames the frames you can use:
"xplot filename.xyz -i ::10 -m"

A list with all the options can be found with 
"xplot -h"
