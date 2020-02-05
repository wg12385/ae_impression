# Makes docker container for running training HPS

FROM continuumio/miniconda3:4.6.14

# Install MKL libraries #######################################################################
RUN apt-get update --yes && apt-get install --yes \
	git \
	make \
	cmake \
	gcc \
	g++ \
	gnupg2 \
	cpio \
	gfortran \
	wget \
	tar \
	libxrender1 \
	libblas-dev \
	liblapack-dev

## If this line fails to find de.debian.org etc, look up local DNS server address then add this to the /etc/docker/daemon.json file like this:
## (137.xxx.x.xx is the server address you found)
#{
#	"dns": ["137.222.0.36", "137.222.0.36"]
#}
# If that still doesnt work, add the following line to /etc/default/docker
# DOCKER_OPTS="--dns 137.222.0.36 --dns 137.222.0.39"
# in either case, restart docker (sudo service docker restart)

# Install MKL #################################################################################

RUN cd /tmp && \
	wget --show-progress --quiet http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/15275/l_mkl_2019.3.199.tgz && \
	tar -xzf l_mkl_2019.3.199.tgz
RUN cd /tmp/l_mkl_2019.3.199 && \
	sed -i 's/ACCEPT_EULA=decline/ACCEPT_EULA=accept/g' silent.cfg && \
	sed -i 's/ARCH_SELECTED=ALL/ARCH_SELECTED=INTEL64/g' silent.cfg && \
	sed -i 's/COMPONENTS=DEFAULTS/COMPONENTS=;intel-comp-l-all-vars__noarch;intel-comp-nomcu-vars__noarch;intel-openmp__x86_64;intel-tbb-libs__x86_64;intel-mkl-common__noarch;intel-mkl-installer-license__noarch;intel-mkl-core__x86_64;intel-mkl-core-rt__x86_64;intel-mkl-doc__noarch;intel-mkl-doc-ps__noarch;intel-mkl-gnu__x86_64;intel-mkl-gnu-rt__x86_64;intel-mkl-common-ps__noarch;intel-mkl-core-ps__x86_64;intel-mkl-common-c__noarch;intel-mkl-core-c__x86_64;intel-mkl-common-c-ps__noarch;intel-mkl-tbb__x86_64;intel-mkl-tbb-rt__x86_64;intel-mkl-gnu-c__x86_64;intel-mkl-common-f__noarch;intel-mkl-core-f__x86_64;intel-mkl-gnu-f-rt__x86_64;intel-mkl-gnu-f__x86_64;intel-mkl-f95-common__noarch;intel-mkl-f__x86_64;intel-mkl-psxe__noarch;intel-psxe-common__noarch;intel-psxe-common-doc__noarch;intel-compxe-pset/g' silent.cfg && \
	./install.sh -s silent.cfg && \
	cd .. && rm -rf * && \
	rm -rf /opt/intel/.*.log /opt/intel/compilers_and_libraries_2019.3.199/licensing && \
	echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/intel.conf && \
	ldconfig && \
	echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> /etc/bash.bashrc

###############################################################################################

# Download aE and qml python libraries ########################################################
RUN git clone https://github.com/wg12385/auto-ENRICH.git opt/autoENRICH
RUN git clone https://github.com/qmlcode/qml.git --branch develop opt/qml && \
	cd opt/qml && git checkout 18f4e9b4fec38b71c9246fb2735a4821229e37b8 && cd ../../

# Create conda environment ####################################################################
RUN conda env create -f opt/autoENRICH/autoENRICH/IMPgen1.yml
SHELL ["conda", "run", "-n", "IMPgen1", "/bin/bash", "-c"]
RUN python -m pip install opt/qml/
RUN python -m pip install opt/autoENRICH/
SHELL ["/bin/bash", "-c"]

RUN echo 'alias autoENRICH="python /opt/autoENRICH/autoENRICH.py" ' >> /.bashrc && \
	echo 'alias IMPRESSION="python /opt/autoENRICH/IMPRESSION.py" ' >> /.bashrc && \
	echo 'alias aEutils="python /opt/autoENRICH/aE_utils.py" ' >> /.bashrc && \
	mkdir output

#ADD RUN RUN
