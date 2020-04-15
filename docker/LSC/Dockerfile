FROM centos:centos7

RUN yum -y install https://download1.rpmfusion.org/free/el/rpmfusion-free-release-7.noarch.rpm
RUN yum -y install epel-release\
           git vim sudo tar unzip wget gcc gcc-c++ make autoconf automake cmake gcc-gfortran\
        python3 python3-pip python3-devel \
        openssh-server initscripts &&\
    yum -y update &&\
    yum clean all



RUN yum -y install libffi-devel \
        freetype6-devel \
        geos \
        geos-devel \
        netcdf-devel \
        libpng12-devel \
        proj\
        proj-devel\

#RUN yum -y install  \
#        openssl-devel \
#        gdal-python \
#        netcdf4-python \
#        python-nose \
#       numpy \
#       python-pyside \
#       sqlite-devel \
#       sqlite\
#       scipy \
#       python-sqlalchemy \
#        python-qt4 \
#        python-yaml \
#        python-matplotlib \
#        python-basemap \
#        texinfo \
#        subversion


RUN pip3 install numexpr

#------------- spatialindex
RUN  cd /tmp  &&\
     wget http://download.osgeo.org/libspatialindex/spatialindex-src-1.8.0.tar.gz &&\
     tar -xvf spatialindex-src-1.8.0.tar.gz &&\
     cd spatialindex-src-1.8.0 &&\
     ./configure --prefix=/usr/local &&\
     make &&\
     make install &&\
     cd /tmp &&\
     rm -rf /tmp/spatialindex-src-1.8.0.tar.gz /tmp/spatialindex-src-1.8.0  &&\
     cp /usr/local/lib/libspatialindex_c.so.3 /usr/lib/ &&\
     ldconfig

#------------- rtree
RUN cd /tmp &&\
    wget https://pypi.python.org/packages/16/68/c79bed9d7fa4479b37d7a898b71d1e89099a8c15fac5c5d94328cb0a1ac9/Rtree-0.8.2.tar.gz &&\
    tar -xvf Rtree-0.8.2.tar.gz &&\
    cd Rtree-0.8.2 &&\
    python3 setup.py install &&\
    cd /tmp &&\
    rm -rf /tmp/Rtree-0.8.2.tar.gz /tmp/Rtree-0.8.2
    

RUN pip3 install pathlib

RUN echo "-----------------Install PySCHISM -----------------" &&\
    cd /tmp &&\
    git clone https://github.com/schism-dev/pyschism &&\
    cd /tmp/pyschism/ &&\
    pip3 install .


RUN echo "-----------------Install PySCHISM -----------------" &&\
    cd /tmp &&\
    git clone https://github.com/metocean/schism-public.git


RUN pip3 install pyyaml numexpr
RUN pip3 install -U pip
RUN yum -y install mesa-libGL
RUN pip3 install PyQt5


RUN cd /usr/local/lib/python3.6/site-packages/pyschism/mesh/ &&\
    sed -i "s|@lru_cache|@lru_cache()|g" base.py &&\
    sed -i "s|@lru_cache|@lru_cache()|g" hgrid.py &&\
    sed -i "s|@lru_cache|@lru_cache()|g" vgrid.py &&\
    sed -i "s|@lru_cache|@lru_cache()|g" mesh.py

RUN cd /usr/local/lib/python3.6/site-packages/pyschism/forcing/ &&\
    sed -i "s|@lru_cache|@lru_cache()|g" tides.py
##
RUN cd /usr/local/lib/python3.6/site-packages/pyschism/driver/&&\
    sed -i "s|@lru_cache|@lru_cache()|g" driver.py &&\
    sed -i "s|@lru_cache|@lru_cache()|g" param.py


RUN yum -y install qt qt-devel libxkbcommon-x11

USER $USER
ENV HOME /home/$USER
    