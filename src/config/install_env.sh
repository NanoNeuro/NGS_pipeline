# With this script we are going to install all the necessary programs as well as the conda environment
# For the necessary libraries.

DIR_PROGRAMS="/home/nanoneuro/Programs"
DIR_DOWNLOADS="/home/nanoneuro/Downloads"


# TODO: INSTALL DOCKERdocker run --rm -it amazon/aws-cli help


# Install apt programs
apt install awscli






# Create conda environment
conda deactivate
conda remove --name NGS_pipeline --all -y
mamba env create --file src/install/conda_env.yaml

conda activate NGS_pipeline







# Install nextflow
cd $DIR_PROGRAMS
wget -qO- https://get.nextflow.io | bash
export PATH="$DIR_PROGRAMS/bin:$PATH"
source ~/.bashrc






# Install centrifuge | we need to install it locally because it fails with the conda installation
DIR_CENTRIFUGE="$DIR_PROGRAMS/centrifuge"

cd $DIR_DOWNLOADS
git clone https://github.com/infphilo/centrifuge
cd centrifuge

sed -i -e 's/“ERROR: please define __cpuid() for this build.\n”/"ERROR: please define __cpuid() for this build.\n"/g' processor_support.h   # Las comillas son raras y dan fallo.
# Si falla, cambiarlo a mano

make
sudo make install prefix=$DIR_CENTRIFUGE

## Add path to bashrc to make the dir recognisable
export PATH="$DIR_CENTRIFUGE/bin:$PATH"
source ~/.bashrc






# Check nextflow is working 
nextflow run nf-core/rnaseq -r 3.13.2 -profile docker,test --outdir results/pruebas/rnaseq -resume
nextflow run nf-core/smrnaseq -r 2.2.4 -profile docker,test --outdir results/pruebas/smrnaseq -resume
nextflow run nf-core/circrna -r dev -profile docker,test --outdir results/pruebas/circrna -resume
nextflow run nf-core/circdna -r 1.0.4 -profile docker,test --outdir results/pruebas/circdna -resume
nextflow run nf-core/taxprofiler -r 1.1.2 -profile docker,test --outdir results/pruebas/taxprofiler -resume
nextflow run nf-core/scrnaseq -r 2.4.1 -profile docker,test --outdir results/pruebas/scrnaseq -resume


