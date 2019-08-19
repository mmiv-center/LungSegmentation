FROM ubuntu:18.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget libboost-filesystem1.62-dev libboost-timer1.62-dev \
    libboost-system1.62-dev

# install itk
RUN cd /tmp/ \
    && wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.0.0/InsightToolkit-5.0.0.tar.gz \
    && cd /opt/ && tar xzvf /tmp/InsightToolkit-5.0.0.tar.gz && cd /opt/InsightToolkit-5.0.0 \
    && mkdir bin && cd bin && cmake .. && make

RUN mkdir /LungSegmentation && cd /LungSegmentation/ \
    && git clone https://github.com/mmiv-center/LungSegmentation.git . \
    && cmake . && make

ENTRYPOINT [ "/LungSegmentation/LungSegmentation" ]
