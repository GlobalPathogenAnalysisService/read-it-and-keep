FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:$PATH
ENV LANG=C.UTF-8

ARG RIK_WF_DIR=/readItAndKeep
RUN mkdir -p $RIK_WF_DIR/.ci/
COPY .ci/install_dependencies.sh $RIK_WF_DIR/.ci/install_dependencies.sh
RUN $RIK_WF_DIR/.ci/install_dependencies.sh /bioinf-tools

COPY . $RIK_WF_DIR
WORKDIR $RIK_WF_DIR
RUN pip3 install tox \
    && cd /readItAndKeep/src \
    && make \
    && make test \
    && cd /bioinf-tools \
    && cp -s /readItAndKeep/src/readItAndKeep .

ENTRYPOINT [ "readItAndKeep" ]