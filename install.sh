path=$(pwd)

if !(grep -q ${path} ~/.bashrc) 
then
echo "export PATH=${path}/render_atoms:\${PATH}" >> ~/.bashrc
echo "export PYTHONPATH=${path}:\${PYTHONPATH}" >> ~/.bashrc
else
echo "render_atoms already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi

chmod +x render_atoms/render_image

source ~/.bashrc
