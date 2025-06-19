path=$(pwd)

if !(grep -q ${path} ~/.bashrc)
then
echo "export PATH=${path}/bin:\${PATH}" >> ~/.bashrc
echo "export PYTHONPATH=${path}:\${PYTHONPATH}" >> ~/.bashrc
else
echo "atomsplot already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi

chmod +x bin/atomsplot

source ~/.bashrc
