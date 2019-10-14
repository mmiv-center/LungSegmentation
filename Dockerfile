FROM ubuntu:18.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget libboost-filesystem1.65-dev libboost-timer1.65-dev \
    libboost-system1.65-dev libboost-date-time1.65-dev

# install itk
RUN cd /tmp/ \
    && wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.0.1/InsightToolkit-5.0.1.tar.gz \
    && cd /opt/ && tar xzvf /tmp/InsightToolkit-5.0.1.tar.gz && cd /opt/InsightToolkit-5.0.1 \
    && mkdir bin && cd bin && cmake .. && make

RUN mkdir /LungSegmentation && cd /LungSegmentation/ \
    && git clone https://github.com/mmiv-center/LungSegmentation.git . \
    && echo "Change this string to make this rebuild on docker build" \
    && cmake . && make

ENTRYPOINT [ "/LungSegmentation/LungSegmentation" ]
