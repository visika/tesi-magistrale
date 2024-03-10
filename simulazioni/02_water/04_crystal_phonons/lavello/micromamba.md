micromamba create -n mace
micromamba activate mace
micromamba install python=3.11 -c conda-forge


# https://pytorch.org/get-started/locally/
# 2.1.2 I hit a bug so use nightly if you want to run on cpu use the cpu version of python
### pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
python3 -m pip install --pre torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cu121
or cuda
# if you use cuda 11.8 use this instead for pythorch
#pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

python3 -m pip install git+https://gitlab.com/ase/ase.git
python3 -m pip install git+https://github.com/ACEsuit/mace.git
python3 -m pip install git+https://github.com/pfnet-research/torch-dftd.git
