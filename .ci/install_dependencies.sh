#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  git \
  liblzma-dev \
  libbz2-dev \
  python3 \
  python3-pip \
  python3-setuptools \
  wget \
  zlib1g-dev

python3 -m pip install pytest pyfastaq

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root


#_________________________ ART ____________________________#
cd $install_root
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
tar xf artbinmountrainier2016.06.05linux64.tgz
rm artbinmountrainier2016.06.05linux64.tgz
cp -s art_bin_MountRainier/art_illumina .


#_________________________ badread ________________________#
cd $install_root
wget -q https://github.com/rrwick/Badread/archive/refs/tags/v0.2.0.tar.gz
tar xf v0.2.0.tar.gz
rm v0.2.0.tar.gz
cd Badread-0.2.0
python3 -m pip install .



