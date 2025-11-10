# Installation of delimtools and software to obtain and analyse testing data

Some test to explain that these instructions will work on Unix.

### Clone delimtools repository, create temp and software directories
##### !!! REQUIRED FOR FULL DELIMTOOLS FUNCTIONALITY !!!
```bash
sudo apt install git
git clone https://github.com/boopsboops/delimtools-testing.git --depth 1
cd delimtools-testing
mkdir temp software

# copy source code for abgd and asap
# original source and backup
# wget https://bioinfo.mnhn.fr/abi/public/asap/last.tgz -O asap.tgz
# wget https://bioinfo.mnhn.fr/abi/public/abgd/last.tgz -O abgd.tgz
cp assets/abgd.tgz software/
cp assets/asap.tgz software/

# wget http://download.r-forge.r-project.org/src/contrib/splits_1.0-20.tar.gz
# wget https://github.com/pedrosenna/drat/blob/master/docs/src/contrib/bGMYC_1.0.3.tar.gz
cp assets/splits_1.0-20.tar.gz software/
cp assets/bGMYC_1.0.3.tar.gz software/

```

### Install R packages to run delimtools-testing
##### !!! REQUIRED FOR FULL DELIMTOOLS FUNCTIONALITY !!!
```bash
# need to be in 'delimtools-testing'
Rscript -e "renv::restore()"

```


### Install mPTP
##### !!! REQUIRED FOR FULL DELIMTOOLS FUNCTIONALITY !!!
```bash
# need to be in 'delimtools-testing/software'
sudo apt-get install libgsl0-dev flex bison autotools-dev autoconf
git clone https://github.com/Pas-Kapli/mptp.git --depth 1
cd mptp
./autogen.sh
./configure
make
cd ..
```


### Install ASAP
##### !!! REQUIRED FOR FULL DELIMTOOLS FUNCTIONALITY !!!
```bash
tar -xzvf asap.tgz
cd ASAP
make
mkdir bin
mv asap bin/
cd ..
```


### Install ABGD
##### !!! REQUIRED FOR FULL DELIMTOOLS FUNCTIONALITY !!!
```bash
tar -xzvf abgd.tgz
cd Abgd
make
mkdir bin
mv abgd bin/
cd ..
```


### Install raxml-ng
##### !!! REQUIRED TO GENERATE DATASET !!!
```bash
wget https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip
unzip raxml-ng_v1.2.2_linux_x86_64.zip -d raxml-ng
echo "export PATH=$(pwd)/raxml-ng:\$PATH" >> ~/.bashrc
exec "$SHELL"
```


### Install BEAST
##### !!! REQUIRED TO GENERATE DATASET !!!
```bash
wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta4/BEAST_X_v10.5.0-beta4.tgz
tar -xzvf BEAST_X_v10.5.0-beta4.tgz
echo "export PATH=$(pwd)/BEASTv10.5.0/bin:\$PATH" >> ~/.bashrc
exec "$SHELL"
cd ..
```


### Install ncbi-supermatrix
##### !!! REQUIRED TO GENERATE DATASET !!!
```bash
# need to be in 'delimtools-testing'
git clone https://github.com/boopsboops/ncbi-supermatrix.git
cd ncbi-supermatrix
# read instructions in README.md to install software required for ncbi-supermatrix
Rscript -e "renv::restore()"
exec "$SHELL"
cd ..
```
